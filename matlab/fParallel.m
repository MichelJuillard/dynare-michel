
function fparallel(fblck,nblck,whoiam,ThisMatlab,fname);

    

    if (ThisMatlab)
        dynareroot = dynare_config();

    end

feval(fname,fblck,nblck,whoiam,ThisMatlab);


%%% Sincronismo "Esterno" %%%%%%%%%%%%%
%%% Ogni Processo quando ha finito lo notifica cancellando un file ... Magari Sistemma 
% keyboard;
if (ThisMatlab)
    if(whoiam)
       load([ fname,'_input'],'MasterName','DyMo','options_' )
        %            fid1 = fopen('P',int2str(whoiam),'End.txt','w+');
        %            fclose(fid1);
        if options_.parallel(ThisMatlab).Local
            
            delete(['P_',fname,'_',int2str(whoiam),'End.txt']);
        else
            delete(['\\',MasterName,'\',DyMo(1),'$\',DyMo(4:end),'\P_',fname,'_',int2str(whoiam),'End.txt']);
%             system (['del \\',MasterName,'\',DyMo(1),'$\',DyMo(4:end),'\P',int2str(whoiam),'End.txt']);
        end
    end

end


if (ThisMatlab==0)
     return;
else
    exit;
end