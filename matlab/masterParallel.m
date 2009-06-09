function [fOutVar,nBlockPerCPU, totCPU] = masterParallel(Parallel,fBlock,nBlock,NamFileInput,fname,fInputVar,fGlobalVar)

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

global options_

totCPU=0;
DyMo=pwd;
fInputVar.DyMo=DyMo;
if isunix,
  [tempo, MasterName]=system(['ifconfig  | grep ''inet addr:''| grep -v ''127.0.0.1'' | cut -d: -f2 | awk ''{ print $1}''']);
else    
  [tempo, MasterName]=system('hostname');
end
MasterName=deblank(MasterName);
fInputVar.MasterName = MasterName;

if exist('fGlobalVar'),
save([fname,'_input.mat'],'fInputVar','fGlobalVar') 
else
save([fname,'_input.mat'],'fInputVar') 
end
save([fname,'_input.mat'],'Parallel','-append') 
for j=1:length(Parallel),
    nCPU(j)=length(Parallel(j).NumCPU);
    totCPU=totCPU+nCPU(j);
end
% keyboard;
nCPU=cumsum(nCPU);
offset0 = fBlock-1;
if (nBlock-offset0)>totCPU,
    diff = mod((nBlock-offset0),totCPU);
    nBlockPerCPU(1:diff) = ceil((nBlock-offset0)/totCPU);
    nBlockPerCPU(diff+1:totCPU) = floor((nBlock-offset0)/totCPU);
else
    nBlockPerCPU(1:nBlock-offset0)=1;
    totCPU = nBlock-offset0;
end

% if totCPU==1,
%     if Parallel.Local == 1
%         State= system (['psexec -W ',DyMo, ' -a ',int2str(Parallel.NumCPU(1)),' -realtime  matlab -nosplash -nodesktop -minimize -r fParallel(',int2str(fBlock),',',int2str(nBlock),',1,1,''',fname,''')']);
%         
%     else
%         if ~strcmp(Parallel.PcName,MasterName),
%             delete(['\\',Parallel(1).PcName,'\',Parallel(1).RemoteDrive,'$\',Parallel(1).RemoteFolder,'\*.*']);
%             adir=ls(['\\',Parallel(1).PcName,'\',Parallel(1).RemoteDrive,'$\',Parallel(1).RemoteFolder,'\']);
%             for j=3:size(adir,1)
%                 rmdir(['\\',Parallel(1).PcName,'\',Parallel(1).RemoteDrive,'$\',Parallel(1).RemoteFolder,'\',adir(j,:)],'s')
%             end
%         end
%         
%         system (['xcopy ',fname,'_input.mat "\\',Parallel.PcName,'\',Parallel.RemoteDrive,'$\',Parallel.RemoteFolder,'" /Y']);
%         for j=1:size(NamFileInput,1)
%             copyfile([NamFileInput{j,1},NamFileInput{j,2}],['\\',Parallel(1).PcName,'\',Parallel(1).RemoteDrive,'$\',Parallel(1).RemoteFolder,'\',NamFileInput{j,1}])
%         end
%         
%               
%         State= system (['psexec \\',Parallel.PcName,' -e -u ',Parallel.user,' -p ',Parallel.passwd,' -W ',Parallel.RemoteDrive,':\',Parallel.RemoteFolder,'\ -a ',int2str(Parallel.NumCPU(1)),' -realtime  matlab -nosplash -nodesktop -minimize -r fParallel(',int2str(fBlock),',',int2str(nBlock),',1,1,''',fname,''')']);
%         
%         system (['xcopy "\\',Parallel.PcName,'\',Parallel.RemoteDrive,'$\',Parallel.RemoteFolder,'\',fname,'_output_1.mat" /Y']);
%         for j=1:size(NamFileOutput,1)
%             copyfile(['\\',Parallel(1).PcName,'\',Parallel(1).RemoteDrive,'$\',Parallel(1).RemoteFolder,'\',NamFileOutput{j,1},NamFileOutput{j,2}],NamFileOutput{j,1})
%         end
%     end
%     load([fname,'_output_1.mat'],'fOutputVar');
%     
% else
    mydelete(['comp_status_',fname,'*.mat'])
    mydelete(['P_',fname,'*End.txt']);
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
          if isunix,
            if exist('OCTAVE_VERSION')
            command1=['octave --eval fParallel\(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',\''',fname,'\''\) &'];
%             command1=['ssh localhost "cd ',DyMo, '; ',matlabroot,'/bin/matlab -nosplash -nodesktop -minimize -r fParallel\(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',\''',fname,'\''\);" &'];
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
%           keyboard;
            if isunix,
              [tempo, RemoteName]=system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' "ifconfig  | grep \''inet addr:\''| grep -v \''127.0.0.1\'' | cut -d: -f2 | awk \''{ print $1}\''"']);
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
                    if isunix,
                      system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' rm -fr ',Parallel(indPC).RemoteFolder,'/*']);
                    else
%                       delete(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\*.*']);
                      mydelete('*.*',['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\']);
                      adir=dir(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\']);
                      for jdir=3:length(adir)
                        STATUS = rmdir(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\',adir(jdir).name],'s');
                        if STATUS == 0,
                          disp(['Warning!: Directory ',adir(jdir).name,' could not be removed from ',Parallel(indPC).PcName,'.'])
                        end
                      end
                    end
                    
                    if isunix,
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
            
            if isunix,
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
    
    if isunix,
    system('sh ConcurrentCommand1.bat &');
    pause(1)
    else
    system('ConcurrentCommand1.bat');
    end
    
    t0=cputime;
    t00=cputime;
    hh=NaN(1,nBlock);
    if exist('OCTAVE_VERSION'),
        diary off;
%         frmt_done = '';
%         for j=1:totCPU,
%           frmt_done = [frmt_done,' %3.f%% '];
%         end
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
%         keyboard;
        %             if (cputime-t0)>10,
        stax = dir(['comp_status_',fname,'*.mat']);
        for j=1:length(stax),
            
            try
                load(stax(j).name)
                pcerdone(j) = prtfrc;
              if exist('OCTAVE_VERSION'),
                statusString = [statusString, waitbarString, ', %3.f%% done! '];
              else
                status_String{j} = waitbarString;  
                status_Title{j} = waitbarTitle;  
                idCPU(j) = njob;
%                 idThisMatlab(j) = ThisMatlab;
%                 idCPU(j) = min((find(cumBlockPerCPU>=njob)));
              end
                if prtfrc==1, delete(stax(j).name), end
            catch
                
            end
%             if ~exist('OCTAVE_VERSION'),
%                 if ishandle(hh(njob)),
%                     waitbar(prtfrc,hh(njob),waitbarString);
%                     if prtfrc==1, close(hh(njob));  end
%                 else
%                     hh(njob) = waitbar(0,waitbarString);
%                     set(hh(njob),'Name',['Parallel ',waitbarTitle]);
%                 end
%             end

        end
            if exist('OCTAVE_VERSION'),
              printf([statusString,'\r'], 100 .* pcerdone);
            else
                  figure(hfigstatus),
              for j=1:length(stax),
                  axes(hstatus(idCPU(j))),
%                   delete(get(gca,'children'))
                  hpat = findobj(hstatus(idCPU(j)),'Type','patch');
                  if ~isempty(hpat),
                    set(hpat,'XData',[0 0 pcerdone(j) pcerdone(j)])
                  else
                    patch([0 0 pcerdone(j) pcerdone(j)],[0 1 1 0],'r','EdgeColor','r')
                  end
%                   xlabel(status_String{j});
                  title([status_Title{j},' - ',status_String{j}]);
              end
            end
        %                 disp(' ')
        %                 t0=cputime;
        %             end
        if isempty(dir(['P_',fname,'_*End.txt'])) 
            mydelete(['comp_status_',fname,'*.mat'])
            if ~exist('OCTAVE_VERSION'),
%                 for j=1:length(hh),
%                     if ishandle(hh(j)),
%                         close(hh(j))
%                     end
%                 end
                close(hfigstatus),
            else
                printf('\n');
                diary on;
            end

%             for j=1:indPC,
%                 if Parallel(j).Local==0 & ~strcmp(Parallel(indPC).PcName,MasterName),
%                     for jfil = 1:size(NamFileOutput,1)
%                     system (['xcopy "\\',Parallel(j).PcName,'\',Parallel(j).RemoteDrive,'$\',Parallel(j).RemoteFolder,'\',NamFileOutput{jfil,1},NamFileOutput{jfil,2},'" ' ,NamFileOutput{jfil,1},' /Y']);
%                     end
%                     system (['xcopy "\\',Parallel(j).PcName,'\',Parallel(j).RemoteDrive,'$\',Parallel(j).RemoteFolder,'\',fname,'_output_*.mat" /Y']);
%                     
%                 end
%             end
            break
        end
    end
    delete([fname,'_input.mat'])
    for j=1:totCPU,
      load([fname,'_output_',int2str(j),'.mat'],'fOutputVar');
      delete([fname,'_output_',int2str(j),'.mat']);
      fOutVar(j)=fOutputVar;

      
    end
      
delete ConcurrentCommand1.bat

% end
