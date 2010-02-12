function [fOutVar,nBlockPerCPU, totCPU] = masterParallel(Parallel,fBlock,nBlock,NamFileInput,fname,fInputVar,fGlobalVar,Parallel_info)
% Top-level function called on the master computer when parallelizing a task.
%
% The number of parallelized threads will be equal to (nBlock-fBlock+1).
%
% INPUTS
%   Parallel [struct vector]   copy of options_.parallel
%   fBlock [int]               index number of the first thread
%                              (between 1 and nBlock)
%   nBlock [int]               index number of the last thread
%   NamFileInput [cell array]  containins the list of input files to be
%                              copied in the working directory of remote slaves
%                              2 columns, as many lines as there are files
%                              - first column contains directory paths
%                              - second column contains filenames
%   fname [string]             name of the function to be parallelized, and
%                              which will be run on the slaves
%   fInputVar [struct]         structure containing local variables to be used
%                              by fName on the slaves
%   fGlobalVar [struct]        structure containing global variables to be used
%                              by fName on the slaves
%
% OUTPUT
%   fOutVar [struct vector]    result of the parallel computation, one
%                              struct per thread
%   nBlockPerCPU [int vector]  for each CPU used, indicates the number of
%                              threads run on that CPU
%   totCPU [int]               total number of CPU used (can be lower than
%                              the number of CPU declared in "Parallel", if
%                              the number of required threads is lower)

% Copyright (C) 2009 Dynare Team
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

% Determine my hostname and my working directory
if ~isempty(Parallel_info)
    if Parallel_info.leaveSlaveOpen,
        [fOutVar,nBlockPerCPU, totCPU] = masterParallelMan(Parallel,fBlock,nBlock,NamFileInput,fname,fInputVar,fGlobalVar);
        return
    end        
end
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

% Save input data for use by the slaves
if exist('fGlobalVar'),
    save([fname,'_input.mat'],'fInputVar','fGlobalVar') 
else
    save([fname,'_input.mat'],'fInputVar') 
end
save([fname,'_input.mat'],'Parallel','-append') 

% Determine the total number of available CPUs, and the number of threads to run on each CPU
[nCPU, totCPU, nBlockPerCPU] = distributeJobs(Parallel, fBlock, nBlock);
offset0 = fBlock-1;

% for j=1:length(Parallel),
%     nCPU(j)=length(Parallel(j).NumCPU);
%     totCPU=totCPU+nCPU(j);
% end
% 
% nCPU=cumsum(nCPU);
% offset0 = fBlock-1;
% if (nBlock-offset0)>totCPU,
%     diff = mod((nBlock-offset0),totCPU);
%     nBlockPerCPU(1:diff) = ceil((nBlock-offset0)/totCPU);
%     nBlockPerCPU(diff+1:totCPU) = floor((nBlock-offset0)/totCPU);
% else
%     nBlockPerCPU(1:nBlock-offset0)=1;
%     totCPU = nBlock-offset0;
% end

% Clean up remnants of previous runs
mydelete(['comp_status_',fname,'*.mat'])
mydelete(['P_',fname,'*End.txt']);

% Create a shell script containing the commands to launch the required tasks on the slaves
fid = fopen('ConcurrentCommand1.bat','w+');
for j=1:totCPU,
    
    indPC=min(find(nCPU>=j));
    
    if indPC>1
        nCPU0 = nCPU(indPC-1);
    else
        nCPU0=0;
    end
    offset = sum(nBlockPerCPU(1:j-1))+offset0;
    
    fid1=fopen(['P_',fname,'_',int2str(j),'End.txt'],'w+');
    fclose(fid1);
    
    if Parallel(indPC).Local == 1, % run on the local machine
        if isunix || (~matlab_ver_less_than('7.4') && ismac),
            if exist('OCTAVE_VERSION')
                command1=['octave --eval fParallel\(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',\''',fname,'\''\) &'];
            else
                command1=['matlab -nosplash -nodesktop -minimize -r fParallel\(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',\''',fname,'\''\) &'];
            end
        else
            if exist('OCTAVE_VERSION')
                command1=['start /B psexec -W ',DyMo, ' -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)),' -low  octave --eval fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
            else
                command1=['start /B psexec -W ',DyMo, ' -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)),' -low  matlab -nosplash -nodesktop -minimize -r fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
            end
        end
    else
        if isunix || (~matlab_ver_less_than('7.4') && ismac),
%             [tempo, RemoteName]=system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' "ifconfig  | grep \''inet addr:\''| grep -v \''127.0.0.1\'' | cut -d: -f2 | awk \''{ print $1}\''"']);
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
                    for jfil=1:size(NamFileInput,1)
                        copyfile([NamFileInput{jfil,1},NamFileInput{jfil,2}],['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\',NamFileInput{jfil,1}])
                    end
                end
            end
        end
        
        if isunix || (~matlab_ver_less_than('7.4') && ismac),
            if exist('OCTAVE_VERSION'),
                command1=['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' "cd ',Parallel(indPC).RemoteFolder, '; octave --eval fParallel\(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',\''',fname,'\''\);" &'];              
            else
                command1=['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' "cd ',Parallel(indPC).RemoteFolder, '; matlab -nosplash -nodesktop -minimize -r fParallel\(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',\''',fname,'\''\);" &'];
            end
        else
            if ~strcmp(Parallel(indPC).PcName,MasterName), % run on a remote machine
                if exist('OCTAVE_VERSION'),
                    command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -u ',Parallel(indPC).user,' -p ',Parallel(indPC).passwd,' -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                              ' -low  octave --eval fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];              
                else
                    command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -u ',Parallel(indPC).user,' -p ',Parallel(indPC).passwd,' -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                              ' -low  matlab -nosplash -nodesktop -minimize -r fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
                end
            else % run on the local machine via the network
                if exist('OCTAVE_VERSION'),
                    command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                              ' -low  octave --eval fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];              
                else
                    command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                              ' -low  matlab -nosplash -nodesktop -minimize -r fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
                end
            end
        end
    end
    fprintf(fid,'%s\n',command1);
end

fclose(fid);

% Run the slaves
if isunix || (~matlab_ver_less_than('7.4') && ismac),
    system('sh ConcurrentCommand1.bat &');
    pause(1)
else
    system('ConcurrentCommand1.bat');
end

% Wait for the slaves to finish their job, and display some progress information meanwhile
t0=cputime;
t00=cputime;
hh=NaN(1,nBlock);
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
            
        end
    end
    if exist('OCTAVE_VERSION'),
        printf([statusString,'\r'], 100 .* pcerdone);
    else
        figure(hfigstatus),
        for j=1:length(stax),
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

% Create return value
for j=1:totCPU,
    load([fname,'_output_',int2str(j),'.mat'],'fOutputVar');
    delete([fname,'_output_',int2str(j),'.mat']);
    fOutVar(j)=fOutputVar;
end

% Cleanup
delete([fname,'_input.mat'])
delete ConcurrentCommand1.bat
