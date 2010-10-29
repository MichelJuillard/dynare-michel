function [ErrorCode] = AnalyseComputationalEnvironment(DataInput, DataInputAdd)
% PARALLEL CONTEXT
% In a parallel context, this function is used to check the user computational request.
% If no error happen the function return 0. The function is able to do it
% for Windows, Linux enviroment for Matlab and Octave software.
% This is a recoursive function. The recursion is on the numbers of
% computer machine in the cluster.
%
% INPUT/OUTPUT description:
%
%
% DataInput
%   is the strcture option_.parallel, with the follow fields:
%
%           Local         Define the computation place: 1 is on local machine, 0 remote
%     ComputerName        Intuitive: contain the computer name.
%           CPUnbr        Intuitive: contain the CPU number.
%         UserName        Intuitive: contain the user name for the ComputerName.
%           Password      Intuitive: contain the password for the user name in ComputerName.
%      RemoteDrive        Drive used for Local/Remote computation (data exchange, etc) must be contain 'RemoteFolder'.
%     RemoteDirectory     Folder in RemoteDrive used for Local/Remote computation.
%     MatlabOctavePath    []
%         DynarePath      []
%
%   This information is typed by the user using the *.mod file,
%   the goal of this function is to check if it correct.
%
%
% DataInputAdd
%   is the structure options_.parallel_info. We use only the string in hte
%   field RemoteTmpFolder (the volatile directory created/destroyed on remote
%   computer). Then for semplicity we derive a string from the struct:

RemoteTmpFolder=DataInputAdd.RemoteTmpFolder;
% DataInputAdd=[];


% The variable ErrorCode is initialized at 0. If there are non problems with
% Local, ComputerName connections,... in general with parallel software execution,
% the ErrorCode is unchanged, in the others cases 1, 2 , ... The values
% table is below.
%
%
%   Table for ErrorCode Values.
%
%   ErrorCode -> 0  Initial Value -> No Error Detected!!!
%   ErrorCode -> > 1  When an error happens. The value 1, 2, 3... are
%   used to specify the kind of error.
%
%   Value 1:    The variable 'Local' has a bad value!
%
%   Value 2:    The variable 'CPUnbr' has a bad value. Parallel Dynare
%               require an input data like [s:d] with s<=d, in this case we
%               have s>d! Or simply the field have no correct length (also
%               empty)!
%         2.1   The user asks to use more CPU of those available.
%         2.2   There are CPU not used!
%
%   Value 3:    The remote computer is unreachable!!!
%
%   Value 4:    The fields user name and/or password are/is empty!
%
%   Value 5:    Remote Drive and/or Remote Folder not exist!
%
%   Value 6:    It is impossible write/read file on remote computer.
%
%   Value 7:    The values user and/or passwd are incorrets or the user have
%               no permissions to execute a Matlab section. Or symply
%               Matlab software is non installed!
%
%   Value 8:    iIt is possible delete remote computational traces!
%
%
%
%
%
%
% Then at the point call of this function it is possible react in a best way, in accord
% with the ErrorCode.

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


ErrorCode=0;

Enviroment=-1;

% Determine a specific operating system or software version when necessary
% for different command (sintax, name, ...).
Enviroment=isunix || (~matlab_ver_less_than('7.4') && ismac);

% Recoursive call:
% ...


disp(' ');
disp(' ');

% The function is composed by two main blocks, determined by the 'Local'
% variable.

% This check can be removed ... in accord with the dynare parser
% strategy.

if ((DataInput.Local == 0) |(DataInput.Local == 1))
    % Continue it is Ok!
    disp('Check on Local Variable ..... Ok!');
    disp(' ');
    disp(' ');
    
else
    disp('The variable "Local" has a bad value!');
    disp(' ');
    disp('ErrorCode 1.');
    disp(' ');
    disp(' ');
    ErrorCode=1;
    return
    
end

%         %%%%%%%%%%  Local (No Network) Computing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Here only the multi-core, or multi-processor avaiable on local
%         machine are involved in parallel computing. No network
%         comunications are required!


% In this case we need to check only the variable 'CPUnbr'.

% We run the parallel code on local computer, so the others fields are automatically
% fixed by Dynare parser. Then the user can also fill them with wrong values.


if (DataInput.Local == 1)
    
    % This check can be removed ... in accord with the dynare parser
    % strategy.
    
    yn=isempty(DataInput.CPUnbr);
    
    if yn==1
        % The field is empty!
        disp('The field "CPUnbr" is empty!');
        disp(' ');
        disp('ErrorCode 2.');
        disp(' ');
        disp(' ');
        ErrorCode=2;
        return
    end
    
    % This check can be removed ... in accord with the dynare parser
    % strategy.
    
    L=length(DataInput.CPUnbr);
    
    if L~=2
        % The field have no correct length!
        disp('The field "CPUnbr" have no length 2!');
        disp(' ');
        disp('ErrorCode 2.');
        disp(' ');
        disp(' ');
        ErrorCode=2;
        return
    end
    
    % This check can be removed ... in accord with the dynare parser
    % strategy.
    
    s=DataInput.CPUnbr(1);
    d=DataInput.CPUnbr(2);
    
    if s>d
        % Bad value s>d!
        disp('In the field "CPUnbr" left side number > right side number!');
        disp(' ');
        disp('ErrorCode 2.');
        disp(' ');
        disp(' ');
        ErrorCode=2;
        return
    end
    
    
    % We look for the information on local computer hardware.
    
    si0=[];
    de0=[];
    
    if Enviroment
        [si0 de0]=system('grep processor /proc/cpuinfo');
    else
        [si0 de0]=system(['psinfo \\']);
    end
    
    
    RealCPUnbr=-1;
    RealCPUnbr=GiveCPUnumber(de0);
    
    % Questo controllo penso che si possa MIGLIORARE!!!!!
    
    if  isempty (RealCPUnbr)
        % An error occurred when we try to know the Cpu/Cores
        % numbers.
        disp('It is impossible determine the number of Cpu/Processor avaiable on this machine!');
        disp(' ');
        disp('ErrorCode 2.');
        disp(' ');
        if Enviroment
            disp('Check the command "$less /proc/cpuinfo" ... !');
        else
            disp('Check if the pstools are installed and are in machine path! And check the command "psinfo \\"');
        end
        disp(' ');
        ErrorCode=2;
        return
    end
    
    
    % Trasforming the input data provided in a form [n1:n2] in a single numerical
    % value.
    
    
    CPUnbrUser=DataInput.CPUnbr(2)-DataInput.CPUnbr(1)+1;
    
    if  CPUnbrUser==RealCPUnbr
        % It is Ok!
        disp('Check on CPUnbr Variable ..... Ok!');
        disp(' ');
        disp(['Hardware have ', num2str(RealCPUnbr),' Cpu/Cores!']);
        disp(['User require ',num2str(CPUnbrUser),' Cpu/Cores!']);
        disp(' ');
        disp(' ');
        
    end
    
    if CPUnbrUser > RealCPUnbr
        disp('Check on CPUnbr Variable ..... Ok!');
        disp(' ');
        disp(['Hardware have ', num2str(RealCPUnbr),' Cpu/Cores!']);
        disp(['User require ',num2str(CPUnbrUser),' Cpu/Cores!']);
        disp(' ');
        disp('Warning! The user asks to use more CPU than those available.');
        disp(' ');
        disp(' ');
        ErrorCode=2.1;
        % return
        
    end
    if CPUnbrUser < RealCPUnbr
        disp('Check on CPUnbr Variable ..... Ok!');
        disp(' ');
        disp(['Hardware have ', num2str(RealCPUnbr),' Cpu/Cores!']);
        disp(['User require ',num2str(CPUnbrUser),' Cpu/Cores!']);
        disp(' ');
        disp('Warning! There are CPU not used!');
        disp(' ');
        disp(' ');
        ErrorCode=2.2;
        % return
    end
    disp('Test for "Local" parallel computation ..... Passed!');
    disp(' ');
    disp(' ');
end



%         %%%%%%%%%%  Cluster Computing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Here we can have many computer with multi-core, or multi-processor avaiable on the
%         network and involved in parallel computing.
%         So in this case we need more sophisticated check.


if (DataInput.Local == 0)
    
    % Now we verify if it is possibile to be connected with the
    % remote computer.
    
    si1=[];
    de1=[];
    
    if Enviroment
        [si1 de1]=system(['ping ', DataInput.ComputerName, ' -c 4']);
    else
        [si1 de1]=system(['ping ', DataInput.ComputerName]);
    end
    
    if (si1)       
        disp(['It is impossibile to be connected to the computer with name "',DataInput.ComputerName,'" using the network!']);
        disp(' ');
        disp('ErrorCode 3.');
        ErrorCode=3;
        disp(' ');
        disp(' ');
        return;
    else
        disp('Check on ComputerName Variable ..... Ok!');
        disp(' ');
        disp(' ');
    end
    
    
    % Now we verify if user name and password are correct and if remote
    % drive and remote folder exist on the remote computer and it is
    % possible to exchange data with them.
    
    if Enviroment
        % This check can be removed ... in accord with the dynare parser
        % strategy.
        
        if (isempty(DataInput.UserName))
            disp('The fields UserName is empty!');
            disp(' ');
            disp('ErrorCode 4.');
            disp(' ');
            disp(' ');
            ErrorCode=4;
            return
        end
        
        % This check can be removed ... in accord with the dynare parser
        % strategy.
        
        if (~isempty(DataInput.Password))
            disp('The fields Password must be empty!');
            disp(' ');
            disp(['Remouve the string ',DataInput.Password,' from this field!']);
            disp(' ');
            disp('ErrorCode 4.');
            disp(' ');
            disp(' ');
            ErrorCode=4;
            return
        end
        disp('Check on UserName Variable ..... Ok!');
        disp(' ');
        disp(' ');
        disp('Check on Password Variable ..... Ok!');
        disp(' ');
        disp(' ');
        
    else
        
        % This check can be removed ... in accord with the dynare parser
        % strategy.
        
        if (isempty(DataInput.UserName)) || (isempty(DataInput.Password))
            disp('The fields UserName and/or Password are/is empty!');
            disp(' ');
            disp('ErrorCode 4.');
            disp(' ');
            disp(' ');
            ErrorCode=4;
            return
        end
        disp('Check on UserName Variable ..... Ok!');
        disp(' ');
        disp(' ');
        disp('Check on Password Variable ..... Ok!');
        disp(' ');
        disp(' ');
        
    end
    
    % Now we very if RemoteDrive and/or RemoteDirectory exist on remote
    % computer!
    
    if Enviroment
        
        % This check can be removed ... in accord with the dynare parser
        % strategy.
        
        if  isempty(DataInput.RemoteDirectory)
            disp('The fields RemoteDirectory is empty!');
            disp(' ');
            disp('ErrorCode 5.');
            disp(' ');
            disp(' ');
            ErrorCode=5;
            return
        end
        
        % This check can be removed ... in accord with the dynare parser
        % strategy.
        
        if (~isempty(DataInput.RemoteDrive))
            disp('The fields RemoteDrive must be empty!');
            disp(' ');
            disp(['Remouve the string ',DataInput.RemoteDrive,' from this field!']);
            disp(' ');
            disp('ErrorCode 5.');
            disp(' ');
            disp(' ');
            ErrorCode=5;
            return
        end
        
        si2=[];
        de2=[];
        
        % Da verificare ... in Linux ...
          [si2 de2]=system(['ssh ',DataInput.UserName,'@',DataInput.ComputerName,' ls ',DataInput.RemoteDirectory,'/',RemoteTmpFolder,'/']);
        
        if (si2)
            disp ('Remote Directory not exist or is not reachables!');
            disp(' ');
            disp('ErrorCode 5.');
            disp(' ');
            disp(' ');
            ErrorCode=5;
            return
        end
        
        disp('Check on RemoteDirectory Variable ..... Ok!');
        disp(' ');
        disp(' ');
        disp('Check on RemoteDrive Variable ..... Ok!');
        disp(' ');
        disp(' ');
        
    else
        % This check can be removed ... in accord with the dynare parser
        % strategy.
        
        if (isempty(DataInput.RemoteDrive)||isempty(DataInput.RemoteDirectory))
            disp('Remote RemoteDrive and/or RemoteDirectory is/are empty!');
            disp(' ');
            disp('ErrorCode 5.');
            disp(' ');
            disp(' ');
            ErrorCode=5;
            return
        end
        
        
        si2=[];
        de2=[];
        [s12 de2]=system(['dir \\',DataInput.ComputerName,'\',DataInput.RemoteDrive,'$\',DataInput.RemoteDirectory,'\',RemoteTmpFolder]);
        
        if (si2)
            disp ('Remote Directory not exist or is not reachables!');
            disp(' ');
            disp('ErrorCode 5.');
            disp(' ');
            disp(' ');
            ErrorCode=5;
            return
        end
        
        disp('Check on RemoteDirectory Variable ..... Ok!');
        disp(' ');
        disp(' ');
        disp('Check on RemoteDrive Variable ..... Ok!');
        disp(' ');
        disp(' ');
        
    end
    
    
    % Now we verify if it possible to exchange data with the remote
    % computer:
    
    % Mettere il percorso in cui si trova 'Tracing.m'!!!
    
    dynareParallelSendFiles('Tracing.m', RemoteTmpFolder,DataInput);
    FindTracing = dynareParallelDir('Tracing.m', RemoteTmpFolder,DataInput);
    
    if (isempty(FindTracing))
        disp ('It is impossible to exchange data with Remote Drive and/or Remote Directory! ErrorCode 6.');
        disp(' ');
        disp('ErrorCode 6.');
        disp(' ');
        disp(' ');
        ErrorCode=6;
        return
    else
        disp('Check on exchange file with remote computer ..... Ok!');
        disp(' ');
        disp(' ');
    end
    
    
    % Now we verify if it is possible execute a matlab/octave section on remote
    % machine when the user is .UserName with password .Password
    
    if exist('OCTAVE_VERSION')
        % OCTAVE!
    else
        % Matlab!
        
        if Enviroment
            % Controllare ... in Linux!
              system(['ssh ',DataInput.UserName,'@',DataInput.ComputerName,' "cd ',DataInput.RemoteDirectory,'/',RemoteTmpFolder,  '; matlab -nosplash -nodesktop -minimize -r Tracing;" &']);
        else
            [NonServeS NenServeD]=system(['start /B /WAIT psexec \\',DataInput.ComputerName,' -e -u ',DataInput.UserName,' -p ',DataInput.Password,' -W ',DataInput.RemoteDrive,':\',DataInput.RemoteDirectory,'\',RemoteTmpFolder ' -low  matlab -nosplash -nodesktop -minimize -r Tracing']);
        end
        
        % Timer da fissare, nei valori di attesa!
        
        t1=fix(clock);
        
        if t1(5)+1>60;
            t2=2;
        else t2=t1(5)+1;
        end
        
        while (1);
            
            nt=fix(clock);
            nt(5)-t2;
            
            if (~isempty (dynareParallelDir('ItIsOk.txt',RemoteTmpFolder,DataInput))) || ((nt(5)-t2)>0)
                if ((nt(5)-t2)>0)
                    ErrorCode=7;
                end
                break
            end
            
        end
        
        if  (ErrorCode==7)
            
            disp ('It is possible execute a matlab section on remote machine!');
            disp(' ');
            disp('ErrorCode 7.');
            disp(' ');
            disp(' ');
            ErrorCode=7;
            return
            
        else
            disp('Check on Matlab execution on remote machine ..... Ok!');
            disp(' ');
            disp(' ');
        end
        
    end
    
    % Now we verify if it is possible delete remote computational traces!
    
    dynareParallelRmDir(RemoteTmpFolder,DataInput);
    
    si3=[];
    de3=[];
    
    if Enviroment
        
        % Da verificare ... in Linux ...
        % [si3 de3]=system(['ssh ',DataInput.UserName,'@',DataInput.ComputerName,' ls ',DataInput.RemoteDirectory,'/']);
    else
        
        [s13 de3]=system(['dir \\',DataInput.ComputerName,'\',DataInput.RemoteDrive,'$\',DataInput.RemoteDirectory,'\',RemoteTmpFolder]);
    end
    
    if isempty(si3)
        disp ('Check on delete remote computational traces ..... Ok!');
        disp(' ');
        disp(' ');
    else
        disp ('It is impossible delete computational traces on remote machine!');
        disp(' ');
        disp('ErrorCode 8.');
        disp(' ');
        disp(' ');
        ErrorCode=8;
        return
    end
    
    
   % Now we check the variable 'CPUnbr'. 
    
    % This check can be removed ... in accord with the dynare parser
    % strategy.
    
    yn=isempty(DataInput.CPUnbr);
    
    if yn==1
        % The field is empty!
        disp('The field "CPUnbr" is empty!');
        disp(' ');
        disp('ErrorCode 2.');
        disp(' ');
        disp(' ');
        ErrorCode=2;
        return
    end
    
    % This check can be removed ... in accord with the dynare parser
    % strategy.
    
    L=length(DataInput.CPUnbr);
    
    if L~=2
        % The field have no correct length!
        disp('The field "CPUnbr" have no length 2!');
        disp(' ');
        disp('ErrorCode 2.');
        disp(' ');
        disp(' ');
        ErrorCode=2;
        return
    end
    
    % This check can be removed ... in accord with the dynare parser
    % strategy.
    
    s=DataInput.CPUnbr(1);
    d=DataInput.CPUnbr(2);
    
    if s>d
        % Bad value s>d!
        disp('In the field "CPUnbr" left side number > right side number!');
        disp(' ');
        disp('ErrorCode 2.');
        disp(' ');
        disp(' ');
        ErrorCode=2;
        return
    end
    
    
    % We look for the information on local computer hardware.
    
    si0=[];
    de0=[];
    
    if Enviroment
        [si0 de0]=system('grep processor /proc/cpuinfo');
    else
        [si0 de0]=system(['psinfo \\']);
    end
    
    
    RealCPUnbr=-1;
    RealCPUnbr=GiveCPUnumber(de0);
    
    % Questo controllo penso che si possa MIGLIORARE!!!!!
    
    if  isempty (RealCPUnbr)
        % An error occurred when we try to know the Cpu/Cores
        % numbers.
        disp('It is impossible determine the number of Cpu/Processor avaiable on this machine!');
        disp(' ');
        disp('ErrorCode 2.');
        disp(' ');
        if Enviroment
            disp('Check the command "$less /proc/cpuinfo" ... !');
        else
            disp('Check if the pstools are installed and are in machine path! And check the command "psinfo \\"');
        end
        disp(' ');
        ErrorCode=2;
        return
    end
    
    
    % Trasforming the input data provided in a form [n1:n2] in a single numerical
    % value.
    
    
    CPUnbrUser=DataInput.CPUnbr(2)-DataInput.CPUnbr(1)+1;
    
    if  CPUnbrUser==RealCPUnbr
        % It is Ok!
        disp('Check on CPUnbr Variable ..... Ok!');
        disp(' ');
        disp(['Hardware have ', num2str(RealCPUnbr),' Cpu/Cores!']);
        disp(['User require ',num2str(CPUnbrUser),' Cpu/Cores!']);
        disp(' ');
        disp(' ');
        
    end
    
    if CPUnbrUser > RealCPUnbr
        disp('Check on CPUnbr Variable ..... Ok!');
        disp(' ');
        disp(['Hardware have ', num2str(RealCPUnbr),' Cpu/Cores!']);
        disp(['User require ',num2str(CPUnbrUser),' Cpu/Cores!']);
        disp(' ');
        disp('Warning! The user asks to use more CPU than those available.');
        disp(' ');
        disp(' ');
        ErrorCode=2.1;
        % return
        
    end
    if CPUnbrUser < RealCPUnbr
        disp('Check on CPUnbr Variable ..... Ok!');
        disp(' ');
        disp(['Hardware have ', num2str(RealCPUnbr),' Cpu/Cores!']);
        disp(['User require ',num2str(CPUnbrUser),' Cpu/Cores!']);
        disp(' ');
        disp('Warning! There are CPU not used!');
        disp(' ');
        disp(' ');
        ErrorCode=2.2;
        % return
    end
    disp('Test for Cluster computation ..... Passed!');
    disp(' ');
    disp(' ');
    
    
    return
end



