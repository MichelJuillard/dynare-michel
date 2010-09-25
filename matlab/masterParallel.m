function [fOutVar,nBlockPerCPU, totCPU] = masterParallel(Parallel,fBlock,nBlock,NamFileInput,fname,fInputVar,fGlobalVar,Parallel_info)
% PARALLEL CONTEXT
% This is the most important function for the management of DYNARE parallel
% computing.
% It is the top-level function called on the master computer when parallelizing a task.


% This function have two main computational startegy for manage the matlab worker (slave process).  
% 0 Simple Close/Open Stategy:
	% In this case the new matlab istances (slave process) are open when
	% necessary and then closed. This can happen many times during the
	% simulation of a model.

% 1 Alway Open Stategy:
	% In this case we have a more sophisticated management of slave processes,
	% which are no longer closed at the end of each job. The slave processes
	% waits for a new job (if exist). If a slave do not receives a new job after a
	% fixed time it is destroyed. This solution removes the computational
	% time necessary to Open/Close new matlab istances.

% The first (point 0) is the default Strategy
% i.e.(Parallel_info.leaveSlaveOpen=0). This value can be changed by the
% user in xxx.mod file or it is changed by the programmer if it necessary to
% reduce the overall computational time. See for example the
% prior_posterior_statistics.m.

% The number of parallelized threads will be equal to (nBlock-fBlock+1).
%
% INPUTS
%  o Parallel [struct vector]   copy of options_.parallel
%  o fBlock [int]               index number of the first thread
%                              (between 1 and nBlock)
%  o nBlock [int]               index number of the last thread
%  o NamFileInput [cell array]  containins the list of input files to be
%                              copied in the working directory of remote slaves
%                              2 columns, as many lines as there are files
%                              - first column contains directory paths
%                              - second column contains filenames
%  o fname [string]             name of the function to be parallelized, and
%                              which will be run on the slaves
%  o fInputVar [struct]         structure containing local variables to be used
%                              by fName on the slaves
%  o fGlobalVar [struct]        structure containing global variables to be used
%                              by fName on the slaves
%
% OUTPUT
%  o fOutVar [struct vector]    result of the parallel computation, one
%                              struct per thread
%  o nBlockPerCPU [int vector]  for each CPU used, indicates the number of
%                              threads run on that CPU
%  o totCPU [int]               total number of CPU used (can be lower than
%                              the number of CPU declared in "Parallel", if
%                              the number of required threads is lower)

% Copyright (C) 2009-2010 Dynare Team
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


% Fix the strategy to be used.
  Strategy=-1;

if ~isempty(Parallel_info)
    Strategy=Parallel_info.leaveSlaveOpen;
end 
for j=1:length(Parallel),
    if isempty(Parallel(j).MatlabPath),
        Parallel(j).MatlabPath = 'matlab';
    end
end
if Strategy==0
   disp('User Strategy Now Is Open/Close (0)');
else
   disp('User Strategy Now Is Always Open (1)');
end
   
if Strategy==1
   % Delete the traces (if exists) of last section computations. 
    persistent initialize
    if isempty(initialize),
       mydelete(['P_slave_*End.txt']);
       mydelete(['slaveParallel_input*.mat']);
       initialize = 0;
       pause(1),
    end
    totCPU=0;        
end

% Determine my hostname and my working directory.
DyMo=pwd;
fInputVar.DyMo=DyMo;
if isunix || (~matlab_ver_less_than('7.4') && ismac) ,
%     [tempo, MasterName]=system(['ifconfig  | grep ''inet addr:''| grep -v ''127.0.0.1'' | cut -d: -f2 | awk ''{ print $1}''']);
    [tempo, MasterName]=system('hostname --fqdn');
else    
    [tempo, MasterName]=system('hostname');
end
MasterName=deblank(MasterName);
fInputVar.MasterName = MasterName;

% Save input data for use by the slaves.
switch Strategy
   case 0
       if exist('fGlobalVar'),
          save([fname,'_input.mat'],'fInputVar','fGlobalVar') 
       else
          save([fname,'_input.mat'],'fInputVar') 
       end
       save([fname,'_input.mat'],'Parallel','-append') 
    
    case 1
       if exist('fGlobalVar'),
          save(['temp_input.mat'],'fInputVar','fGlobalVar')
       else
          save(['temp_input.mat'],'fInputVar')
       end
       save(['temp_input.mat'],'Parallel','-append')
end


% Determine the total number of available CPUs, and the number of threads
% to run on each CPU.

   [nCPU, totCPU, nBlockPerCPU] = distributeJobs(Parallel, fBlock, nBlock);
   offset0 = fBlock-1;  
  

% Clean up remnants of previous runs.
mydelete(['comp_status_',fname,'*.mat'])
mydelete(['P_',fname,'*End.txt']);

% Create a shell script containing the commands to launch the required tasks on the slaves
fid = fopen('ConcurrentCommand1.bat','w+');
for j=1:totCPU,
    
    
    if Strategy==1
        command1 = ' ';
    end
   
    
    indPC=min(find(nCPU>=j));
    
    if indPC>1
        nCPU0 = nCPU(indPC-1);
    else
        nCPU0=0;
    end
    offset = sum(nBlockPerCPU(1:j-1))+offset0;
    
    fid1=fopen(['P_',fname,'_',int2str(j),'End.txt'],'w+');
    fclose(fid1);
    
    if Strategy==1
       fblck = offset+1;
       nblck = sum(nBlockPerCPU(1:j));
       save temp_input fblck nblck fname -append;
       copyfile('temp_input.mat',['slaveJob',int2str(j),'.mat'])
       if Parallel(indPC).Local ==0,
            fid1=fopen(['stayalive',int2str(j),'.txt'],'w+');
            fclose(fid1);
            if isunix || (~matlab_ver_less_than('7.4') && ismac),
                system(['scp stayalive',int2str(j),'.txt ',Parallel(indPC).user,'@',Parallel(indPC).PcName,':',Parallel(indPC).RemoteFolder]);
            else
            copyfile(['stayalive',int2str(j),'.txt'], ['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder]);
            end
            mydelete(['stayalive',int2str(j),'.txt'],'w+');
       end
       % Wait for possibly local alive CPU to start the new job or close by
       % internal criteria.
        pause(1); 
        newInstance = 0;
        
        % Check if j CPU is already alive.
        if isempty( dir(['P_slave_',int2str(j),'End.txt'])); 
            fid1=fopen(['P_slave_',int2str(j),'End.txt'],'w+');
            fclose(fid1);
            newInstance = 1;
            storeGlobalVars( ['slaveParallel_input',int2str(j)]);
            save( ['slaveParallel_input',int2str(j)],'Parallel','-append');
            % Prepare global vars for Slave.
        end
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % The following 'switch - case' code is the core of this function!
   switch Strategy
    case 0
      
        if Parallel(indPC).Local == 1, %Run on the local machine (localhost).
            if isunix || (~matlab_ver_less_than('7.4') && ismac),
                if exist('OCTAVE_VERSION')
                    command1=['octave --eval "addpath(''',Parallel(indPC).DynarePath,'''), fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')" &'];
                else
                    command1=[Parallel(indPC).MatlabPath,' -nosplash -nodesktop -minimize -r "addpath(''',Parallel(indPC).DynarePath,'''), fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')" &'];
                end
            else
                if exist('OCTAVE_VERSION')
                    command1=['start /B psexec -W ',DyMo, ' -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)),' -low  octave --eval "addpath(''',Parallel(indPC).DynarePath,'''), fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')"'];
                else
                    command1=['start /B psexec -W ',DyMo, ' -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)),' -low  ',Parallel(indPC).MatlabPath,' -nosplash -nodesktop -minimize -r "addpath(''',Parallel(indPC).DynarePath,'''), fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')"'];
                end
            end
        else % Parallel(indPC).Local==0: Run using network on remote machine or also on local machine.
            if isunix || (~matlab_ver_less_than('7.4') && ismac),
                %[tempo, RemoteName]=system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' "ifconfig  | grep \''inet addr:\''| grep -v \''127.0.0.1\'' | cut -d: -f2 | awk \''{ print $1}\''"']);
                [tempo, RemoteName]=system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' "hostname --fqdn"']);
                RemoteName=RemoteName(1:end-1);
                RemoteFolder = Parallel(indPC).RemoteFolder;
            else    
                RemoteName = Parallel(indPC).PcName;
                RemoteFolder = [Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder];
            end
        
            remoteFlag=1;

            if strcmpi(RemoteName,MasterName),
                if ~copyfile(['P_',fname,'_',int2str(j),'End.txt'],RemoteFolder),
                    remoteFlag=0;
                end
            end
            if remoteFlag,
                if j==nCPU0+1,
                    if isunix || (~matlab_ver_less_than('7.4') && ismac),
                       system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' rm -fr ',Parallel(indPC).RemoteFolder,'/*']);
                    else
                        mydelete('*.*',['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\']);
                        adir=dir(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\']);
                        for jdir=3:length(adir)
                            STATUS = rmdir(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\',adir(jdir).name],'s');
                            if STATUS == 0,
                                disp(['Warning!: Directory ',adir(jdir).name,' could not be removed from ',Parallel(indPC).PcName,'.'])
                            end
                        end
                    end
              if isunix || (~matlab_ver_less_than('7.4') && ismac),
                    system(['scp ',fname,'_input.mat ',Parallel(indPC).user,'@',Parallel(indPC).PcName,':',Parallel(indPC).RemoteFolder]);
                    for jfil=1:size(NamFileInput,1)
                        if ~isempty(NamFileInput{jfil,1})
                            system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' mkdir -p ',Parallel(indPC).RemoteFolder,'/',NamFileInput{jfil,1}])
                        end
                        system(['scp ',NamFileInput{jfil,1},NamFileInput{jfil,2},' ',Parallel(indPC).user,'@',Parallel(indPC).PcName,':',Parallel(indPC).RemoteFolder,'/',NamFileInput{jfil,1}]);
                    end
              else
                  copyfile([fname,'_input.mat'], ['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder]);
                  for jfil=1:size(NamFileInput,1),
                      if ~isempty(NamFileInput{jfil,1})
                         mkdir(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\',NamFileInput{jfil,1}]);
                      end
                      copyfile([NamFileInput{jfil,1},NamFileInput{jfil,2}],['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\',NamFileInput{jfil,1}])
                  end
              end
            end
        end
        
        if isunix || (~matlab_ver_less_than('7.4') && ismac),
            if exist('OCTAVE_VERSION'),
                command1=['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' "cd ',Parallel(indPC).RemoteFolder, '; octave --eval \"addpath(''',Parallel(indPC).DynarePath,'''), fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''');\" " &'];              
            else
                command1=['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' "cd ',Parallel(indPC).RemoteFolder, '; ',Parallel(indPC).MatlabPath,' -nosplash -nodesktop -minimize -r \"addpath(''',Parallel(indPC).DynarePath,'''), fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''');\" " &'];
            end
        else
            if ~strcmp(Parallel(indPC).PcName,MasterName), % Run on a remote machine!
                if exist('OCTAVE_VERSION'),
                    command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -u ',Parallel(indPC).user,' -p ',Parallel(indPC).passwd,' -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                              ' -low  octave --eval "addpath(''',Parallel(indPC).DynarePath,'''), fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')"'];              
                else
                    command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -u ',Parallel(indPC).user,' -p ',Parallel(indPC).passwd,' -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                              ' -low  ',Parallel(indPC).MatlabPath,' -nosplash -nodesktop -minimize -r "addpath(''',Parallel(indPC).DynarePath,'''), fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')"'];
                end
            else % Run on the local machine via the network
                if exist('OCTAVE_VERSION'),
                    command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                              ' -low  octave --eval "addpath(''',Parallel(indPC).DynarePath,'''), fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')"'];              
                else
                    command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                              ' -low  ',Parallel(indPC).MatlabPath,' -nosplash -nodesktop -minimize -r "addpath(''',Parallel(indPC).DynarePath,'''), fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')"'];
                end
            end
        end
        end
   
    
   case 1
        if Parallel(indPC).Local == 1 & newInstance, % run on the local machine
            if isunix || (~matlab_ver_less_than('7.4') && ismac),
                if exist('OCTAVE_VERSION')
                   %command1=['octave --eval fParallel\(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',\''',fname,'\''\) &'];
                    command1=['octave --eval "addpath(''',Parallel(indPC).DynarePath,'''), slaveParallel(',int2str(j),',',int2str(indPC),')" &'];
                else
                    %command1=[Parallel(indPC).MatlabPath,' -nosplash -nodesktop -minimize -r fParallel\(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',\''',fname,'\''\) &'];
                    command1=[Parallel(indPC).MatlabPath,' -nosplash -nodesktop -minimize -r "addpath(''',Parallel(indPC).DynarePath,'''), slaveParallel(',int2str(j),',',int2str(indPC),')" &'];
                end
            else
                if exist('OCTAVE_VERSION')
                    %command1=['start /B psexec -W ',DyMo, ' -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)),' -low  octave --eval fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
                    command1=['start /B psexec -W ',DyMo, ' -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)),' -low  octave --eval "addpath(''',Parallel(indPC).DynarePath,'''), slaveParallel(',int2str(j),',',int2str(indPC),')"'];
                else
                    %command1=['start /B psexec -W ',DyMo, ' -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)),' -low  ',Parallel(indPC).MatlabPath,' -nosplash -nodesktop -minimize -r fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
                    command1=['start /B psexec -W ',DyMo, ' -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)),' -low  ',Parallel(indPC).MatlabPath,' -nosplash -nodesktop -minimize -r "addpath(''',Parallel(indPC).DynarePath,'''), slaveParallel(',int2str(j),',',int2str(indPC),')"'];
                end
            end
            
        elseif Parallel(indPC).Local==0 %Run using network on remote machine or also on local machine.
        if isunix || (~matlab_ver_less_than('7.4') && ismac),
            %[tempo, RemoteName]=system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' "ifconfig  | grep \''inet addr:\''| grep -v \''127.0.0.1\'' | cut -d: -f2 | awk \''{ print $1}\''"']);
            [tempo, RemoteName]=system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' "hostname --fqdn"']);
            RemoteName=RemoteName(1:end-1);
            RemoteFolder = Parallel(indPC).RemoteFolder;
        else
            RemoteName = Parallel(indPC).PcName;
            RemoteFolder = [Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder];
        end
        
        remoteFlag=1;
        
        if strcmpi(RemoteName,MasterName),
            if ~copyfile(['P_',fname,'_',int2str(j),'End.txt'],RemoteFolder),
                remoteFlag=0;
            end
        end
        if remoteFlag,
            if (j==nCPU0+1) & newInstance, % Clean remote folder!
                if isunix || (~matlab_ver_less_than('7.4') && ismac),
                    system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' rm -fr ',Parallel(indPC).RemoteFolder,'/*']);
                else
                    mydelete('*.*',['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\']);
                    adir=dir(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\']);
                    for jdir=3:length(adir)
                        STATUS = rmdir(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\',adir(jdir).name],'s');
                        if STATUS == 0,
                           disp(['Warning!: Directory ',adir(jdir).name,' could not be removed from ',Parallel(indPC).PcName,'.'])
                        end
                    end
                end
            end
            
            if isunix || (~matlab_ver_less_than('7.4') && ismac),
                % system(['scp ',fname,'_input.mat ',Parallel(indPC).user,'@',Parallel(indPC).PcName,':',Parallel(indPC).RemoteFolder]);
                for jfil=1:size(NamFileInput,1)
                    if ~isempty(NamFileInput{jfil,1})
                        system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' mkdir -p ',Parallel(indPC).RemoteFolder,'/',NamFileInput{jfil,1}])
                    end
                    system(['scp ',NamFileInput{jfil,1},NamFileInput{jfil,2},' ',Parallel(indPC).user,'@',Parallel(indPC).PcName,':',Parallel(indPC).RemoteFolder,'/',NamFileInput{jfil,1}]);
                end
                system(['scp slaveJob',int2str(j),'.mat ',Parallel(indPC).user,'@',Parallel(indPC).PcName,':',Parallel(indPC).RemoteFolder]);
                if newInstance,
                    system(['scp slaveParallel_input',int2str(j),'.mat ',Parallel(indPC).user,'@',Parallel(indPC).PcName,':',Parallel(indPC).RemoteFolder]);
                end
            else
                %copyfile([fname,'_input.mat'], ['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder])
                for jfil=1:size(NamFileInput,1)
                    if ~isempty(NamFileInput{jfil,1})
                        mkdir(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\',NamFileInput{jfil,1}]);
                    end
                    copyfile([NamFileInput{jfil,1},NamFileInput{jfil,2}],['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\',NamFileInput{jfil,1}])
                end
                copyfile(['slaveJob',int2str(j),'.mat'], ['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder]);
                if newInstance,
                    copyfile(['slaveParallel_input',int2str(j),'.mat'], ['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder]);
                end
            end
        end
        
        if newInstance,
            if isunix || (~matlab_ver_less_than('7.4') && ismac),
                if exist('OCTAVE_VERSION'),
                    command1=['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' "cd ',Parallel(indPC).RemoteFolder, '; octave --eval \"addpath(''',Parallel(indPC).DynarePath,'''), slaveParallel(',int2str(j),',',int2str(indPC),');\" " &'];
                else
                    command1=['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' "cd ',Parallel(indPC).RemoteFolder, '; ',Parallel(indPC).MatlabPath,' -nosplash -nodesktop -minimize -r \"addpath(''',Parallel(indPC).DynarePath,'''), slaveParallel(',int2str(j),',',int2str(indPC),');\" " &'];
                end
            else
                if ~strcmp(Parallel(indPC).PcName,MasterName), % Run on a remote machine.
                    if exist('OCTAVE_VERSION'),
                        command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -u ',Parallel(indPC).user,' -p ',Parallel(indPC).passwd,' -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                            ' -low  octave --eval "addpath(''',Parallel(indPC).DynarePath,'''), slaveParallel(',int2str(j),',',int2str(indPC),')"'];
                    else
                        command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -u ',Parallel(indPC).user,' -p ',Parallel(indPC).passwd,' -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                            ' -low  ',Parallel(indPC).MatlabPath,' -nosplash -nodesktop -minimize -r "addpath(''',Parallel(indPC).DynarePath,'''), slaveParallel(',int2str(j),',',int2str(indPC),')"'];
                    end
                else % Run on the local machine via the network.
                    if exist('OCTAVE_VERSION'),
                        command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                            ' -low  octave --eval "addpath(''',Parallel(indPC).DynarePath,'''), slaveParallel(',int2str(j),',',int2str(indPC),')"'];
                    else
                        command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                            ' -low  ',Parallel(indPC).MatlabPath,' -nosplash -nodesktop -minimize -r "addpath(''',Parallel(indPC).DynarePath,'''), slaveParallel(',int2str(j),',',int2str(indPC),')"'];
                    end
                end
            end
        end
    end  
   end

 fprintf(fid,'%s\n',command1);


end
 fclose(fid);

% Run the slaves.
if isunix || (~matlab_ver_less_than('7.4') && ismac),
    system('sh ConcurrentCommand1.bat &');
    pause(1)
else
    system('ConcurrentCommand1.bat');
end

% Wait for the slaves to finish their job, and display some progress
% information meanwhile.


if exist('OCTAVE_VERSION'),
    diary off;
    printf('\n');
else
    hfigstatus = figure('name',['Parallel ',fname],...
                        'MenuBar', 'none', ...
                        'NumberTitle','off');
    vspace = 0.1;
    ncol = ceil(totCPU/10);
    hspace = 0.9/ncol;
    for j=1:totCPU,
        jrow = mod(j-1,10)+1;
        jcol = ceil(j/10);  
        hstatus(j) = axes('position',[0.05/ncol+(jcol-1)/ncol 0.92-vspace*(jrow-1) 0.9/ncol 0.03], ...
                          'box','on','xtick',[],'ytick',[],'xlim',[0 1],'ylim',[0 1]);
    end
    cumBlockPerCPU = cumsum(nBlockPerCPU);
end
pcerdone = NaN(1,totCPU);
while (1)
    
    waitbarString = '';
    statusString = '';
    pause(1)
    
    stax = dir(['comp_status_',fname,'*.mat']);
    
    for j=1:length(stax),   
        try
            load(stax(j).name)
            pcerdone(j) = prtfrc;
            if exist('OCTAVE_VERSION'),
                statusString = [statusString, int2str(j), ' %3.f%% done! '];
            else
                status_String{j} = waitbarString;  
                status_Title{j} = waitbarTitle;  
                idCPU(j) = njob;
            end
            if prtfrc==1, delete(stax(j).name), end
        catch
            if j>1
               j=j-1
            end
        end
    end
    if exist('OCTAVE_VERSION'),
        printf([statusString,'\r'], 100 .* pcerdone);
    else
        figure(hfigstatus),
            for j=1:length(stax)
            try
            axes(hstatus(idCPU(j))),
            hpat = findobj(hstatus(idCPU(j)),'Type','patch');
            if ~isempty(hpat),
                set(hpat,'XData',[0 0 pcerdone(j) pcerdone(j)])
            else
                patch([0 0 pcerdone(j) pcerdone(j)],[0 1 1 0],'r','EdgeColor','r')
            end
            title([status_Title{j},' - ',status_String{j}]);
            catch
                    if j>1
                        j=j-1
                    end
            end
        end
    end
    if isempty(dir(['P_',fname,'_*End.txt'])) 
        mydelete(['comp_status_',fname,'*.mat'])
        if ~exist('OCTAVE_VERSION'),
            close(hfigstatus),
        else
            printf('\n');
            diary on;
        end
        break
    end
end

% Create return value.
for j=1:totCPU,
    load([fname,'_output_',int2str(j),'.mat'],'fOutputVar');
    delete([fname,'_output_',int2str(j),'.mat']);
    fOutVar(j)=fOutputVar;
end

% Cleanup.
switch Strategy
   case 0
       delete([fname,'_input.mat'])
   case 1
       delete(['temp_input.mat'])
end
        
delete ConcurrentCommand1.bat
