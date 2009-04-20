function [nBlockPerCPU, totCPU] = masterParallel(Parallel,nBlock,NamFileInput,NamFileOutput,fname)
global options_

totCPU=0;
DyMo=pwd;
[tempo, MasterName]=dos('hostname');
MasterName=MasterName(1:end-1);
save([fname,'_input.mat'],'DyMo','MasterName','-append') 

for j=1:length(Parallel),
    nCPU(j)=length(Parallel(j).NumCPU);
    totCPU=totCPU+nCPU(j);
end
nCPU=cumsum(nCPU);
if nBlock>totCPU,
    diff = mod(nBlock,totCPU);
    nBlockPerCPU(1:diff) = ceil(nBlock/totCPU);
    nBlockPerCPU(diff+1:totCPU) = floor(nBlock/totCPU);
else
    nBlockPerCPU(1:nBlock)=1;
    totCPU = nBlock;
end

if totCPU==1,
    if Parallel.Local == 1
        State= system (['psexec -W ',DyMo, ' -a ',int2str(Parallel.NumCPU(1)),' -realtime  matlab -nosplash -nodesktop -minimize -r fParallel(1,',int2str(nBlock),',1,1,''',fname,''')']);
        
    else
        if ~strcmp(Parallel.PcName,MasterName),
            delete(['\\',Parallel(1).PcName,'\',Parallel(1).RemoteDrive,'$\',Parallel(1).RemoteFolder,'\*.*']);
            adir=ls(['\\',Parallel(1).PcName,'\',Parallel(1).RemoteDrive,'$\',Parallel(1).RemoteFolder,'\']);
            for j=3:size(adir,1)
                rmdir(['\\',Parallel(1).PcName,'\',Parallel(1).RemoteDrive,'$\',Parallel(1).RemoteFolder,'\',adir(j,:)],'s')
            end
        end
        
        system (['xcopy ',fname,'_input.mat "\\',Parallel.PcName,'\',Parallel.RemoteDrive,'$\',Parallel.RemoteFolder,'" /Y']);
        for j=1:size(NamFileInput,1)
            copyfile([NamFileInput{j,1},NamFileInput{j,2}],['\\',Parallel(1).PcName,'\',Parallel(1).RemoteDrive,'$\',Parallel(1).RemoteFolder,'\',NamFileInput{j,1}])
        end
        
              
        State= system (['psexec \\',Parallel.PcName,' -e -u ',Parallel.user,' -p ',Parallel.passwd,' -W ',Parallel.RemoteDrive,':\',Parallel.RemoteFolder,'\ -a ',int2str(Parallel.NumCPU(1)),' -realtime  matlab -nosplash -nodesktop -minimize -r fParallel(1,',int2str(nBlock),',1,1,''',fname,''')']);
        
        system (['xcopy "\\',Parallel.PcName,'\',Parallel.RemoteDrive,'$\',Parallel.RemoteFolder,'\',fname,'_output_1.mat" /Y']);
        for j=1:size(NamFileOutput,1)
            copyfile(['\\',Parallel(1).PcName,'\',Parallel(1).RemoteDrive,'$\',Parallel(1).RemoteFolder,'\',NamFileOutput{j,1},NamFileOutput{j,2}],NamFileOutput{j,1})
        end
    end
    out=load([fname,'_output.mat']);
    
else
    delete(['comp_status_',fname,'*.mat'])
    fid = fopen('ConcurrentCommand1.bat','w+');
    for j=1:totCPU,
        
        indPC=min(find(nCPU>=j));
        
        if indPC>1
            nCPU0 = nCPU(indPC-1);
        else
            nCPU0=0;
        end
        offset = sum(nBlockPerCPU(1:j-1));
        
        fid1=fopen(['P_',fname,'_',int2str(j),'End.txt'],'w+');
        fclose(fid1);
        
        if Parallel(indPC).Local == 1
%             command1=['start /B /affinity ',int2str(Parallel(indPC).NumCPU(j-nCPU0)+1),' /normal  matlab -nosplash -nodesktop -minimize -sd ',DyMo, ' -r fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
            command1=['start /B psexec -W ',DyMo, ' -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)),' -high  matlab -nosplash -nodesktop -minimize -r fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
            
        else
            
            if ~strcmp(Parallel(indPC).PcName,MasterName),
                if j==nCPU0+1,
                    delete(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\*.*']);
                    adir=ls(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\']);
                    for jdir=3:size(adir,1)
                        rmdir(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\',adir(jdir,:)],'s')
                    end
                    system (['xcopy ',fname,'_input.mat "\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'" /Y']);
                    for jfil=1:size(NamFileInput,1)
                       copyfile([NamFileInput{jfil,1},NamFileInput{jfil,2}],['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\',NamFileInput{jfil,1}])
                    end
                end
              command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -u ',Parallel(indPC).user,' -p ',Parallel(indPC).passwd,' -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                  ' -normal  matlab -r fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
%                   ' -normal  matlab -nosplash -nodesktop -minimize -r fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
            else
              command1=['start /B psexec \\',Parallel(indPC).PcName,' -e -W ',Parallel(indPC).RemoteDrive,':\',Parallel(indPC).RemoteFolder,'\ -a ',int2str(Parallel(indPC).NumCPU(j-nCPU0)), ...
                  ' -normal  matlab -nosplash -nodesktop -minimize -r fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),',''',fname,''')'];
            end
            
            
            
            %                 woDy=pwd;
            
            
            %                 State= system (['psexec \\',Parallel.PcName,' -i -W ',Parallel.RemoteDrive,':\',Parallel.RemoteFolder,'\ -a 0 -realtime  matlab -nosplash -nodesktop -minimize -r fParallel(1,',int2str(nBlock),',0,1)']);
            %                     ' -normal  matlab -nosplash -nodesktop -minimize -r fParallel(',int2str(offset+1),',',int2str(sum(nBlockPerCPU(1:j))),',',int2str(j),',',int2str(indPC),')'];
            %                 command2=['start /MIN psexec \\',Parallel.PcName,' -W ',  DyMo,  ' -a 1 -high matlab -nosplash -nodesktop -minimize -r fParallel(2,2,2,1)'];
            %                 command3=['start /MIN psexec \\',Parallel.PcName,' -W ',  DyMo,  ' -a 2 -high matlab -nosplash -nodesktop -minimize -r fParallel(3,3,3,1)'];
            %                 command4=['start /wait /MIN psexec \\',Parallel.PcName,' -W ', DyMo, ' -a 3 -high matlab -nosplash -nodesktop -minimize -r fParallel(4,4,4,1)'];
            
            %                 State= system (['psexec \\',Parallel.PcName,' -i -W ',Parallel.RemoteDrive,':\',Parallel.RemoteFolder,'\ -a 0 -realtime  matlab -nosplash -nodesktop -minimize -r fParallel(1,',int2str(nBlock),',0,1)']);
            
        end
        fprintf(fid,'%s\n',command1);
    end
    
    fclose(fid);
    
    dos('ConcurrentCommand1.bat')
    delete ConcurrentCommand1.bat
    
    status_string={'Starting remote parrallel computation ... '};
    t0=cputime;
    t00=cputime;
    hh=NaN(1,nBlock);
    while (1)
        pause(1)
        %             if (cputime-t0)>10,
        stax = ls(['comp_status_',fname,'*.mat']);
        for j=1:size(stax,1),
            
            try
                load(stax(j,:))
                %                         status_string{j}=(['Chain ',int2str(whoiam),' at ',num2str(100*jstatus/nruns(whoiam)),'% accept. rate ',num2str(isux/jstatus,4),'.']);
            catch
                
            end
            %                   disp(status_string{j})
            prtfrc = jstatus/nruns(b);
            if ishandle(hh(b)),
                waitbar(prtfrc,hh(b),[ '(' int2str(b) '/' int2str(options_.mh_nblck) ') ' sprintf('%f done, acceptation rate %f',prtfrc,isux/jstatus)]);
                if prtfrc==1, close(hh(b)); delete(stax(j,:)), end
            else
                hh(b) = waitbar(0,['Please wait... Metropolis-Hastings (' int2str(b) '/' int2str(options_.mh_nblck) ')...']);
                set(hh(b),'Name',['Remote Metropolis-Hastings']);
            end
            
        end
        %                 disp(' ')
        %                 t0=cputime;
        %             end
        if isempty(ls(['P_',fname,'_*End.txt'])) 
            delete(['comp_status_',fname,'*.mat'])
            for j=1:length(hh),
                if ishandle(hh(j)),
                    close(hh(j))
                end
            end

            for j=1:indPC,
                if Parallel(j).Local==0 & ~strcmp(Parallel(indPC).PcName,MasterName),
                    for jfil = 1:size(NamFileOutput,1)
                    system (['xcopy "\\',Parallel(j).PcName,'\',Parallel(j).RemoteDrive,'$\',Parallel(j).RemoteFolder,'\',NamFileOutput{jfil,1},NamFileOutput{jfil,2},'" ' ,NamFileOutput{jfil,1},' /Y']);
                    end
                    system (['xcopy "\\',Parallel(j).PcName,'\',Parallel(j).RemoteDrive,'$\',Parallel(j).RemoteFolder,'\',fname,'_output_*.mat" /Y']);
                    
                end
            end
            break
        end
    end
    
end
