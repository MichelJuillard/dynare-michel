function [ErrorCode] = AnalyseComputationalEnviroment(DataInput)
% Input/Output description:
%
% DataInput is the strcture option_.parallel, with the follow fields:
%
%           Local         Define the computation place: 1 is on local machine, 0 remote
%           PcName        Intuitive: contain the computer name.
%           NumCPU        Intuitive: contain the CPU number.
%             user        Intuitive: contain the use name for the PcName.
%           passwd        Intuitive: contain the password for the user name in PcName.
%      RemoteDrive        Drive used for Local/Remote computation (data exchange, etc) must be contain 'RemoteFolder'.
%     RemoteFolder        Folder in RemoteDrive used for Local/Remote computation.
%
%   This information is typed by the user using the *.mod file, 
%   the goal of this function is to check if it correct.
%
%
% The variable ErrorCode is initialized at 0. If there are non problems with 
% Local, PcName connections,... in general with parallel software execution, 
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
%   Value 2:    The variable 'NumCPU' has a bad value. Parallel Dynare
%               require an input data like [s:d] with s<=d, in this case we
%               have s>d!
%         2.1   The user asks to use more CPU of those available.
%         2.2   There are CPU not used!
%
%   Value 3:    The remote computer is unreachable!!!
%
%   Value 4:    The user name and/or password is/are incorrect on the
%               remote computer!
%
%   Value 5:    It is impossible write/read file on remote computer.
%
% Then at the point call of this function it is possible react in a best way, in accord
% with the ErrorCode.

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

ErrorCode=0;


% The function is composed by two main blocks, determined by the 'Local'
% variable. 

if ((DataInput.Local == 0) |(DataInput.Local == 1))
    % Continue it is Ok!
else
    ErrorCode=1;
    return

end


%%%%%%%%%%  Local Machine   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this case we need to check only the variable 'NumCPU'. 

% We run the parallel code on local computer, so the others fields are automatically 
% fixed by Dynare. Then the user can also fill them with wrong values.

if (DataInput.Local == 1)
    
    yn=isempty(DataInput.NumCPU);
    
    if yn==1
        ErrorCode=2;
        return
    end
    
    % We look for the information on local computer hardware.
    
    si=[];
    de=[];
    
    [si de]=system(['psinfo \\']);
    
    RealNumCPU=-1;
    RealNumCPU=GiveCPUnumber(de);
    
    % Trasforming the input data provided in a form [n1:n2] in a single numerical
    % value.   
    
    DataInput.NumCPU=length(DataInput.NumCPU);
    
    if DataInput.NumCPU  == RealNumCPU
        % It is Ok!
    end
    
    if DataInput.NumCPU > RealNumCPU
        ErrorCode=2.1;
        
    end
    if DataInput.NumCPU < RealNumCPU
        ErrorCode=2.2;
    end    
end   

%%%%%%%%%%  Remote Machine   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this case we need more sophisticated check. 


if (DataInput.Local == 0)
    
    si=[];
    de=[]; 
    
    [si de]=system(['ping ', DataInput.PcName]);
    
    if si==1 
        % It is impossiblie to be connected to the
        % remote computer.
        
        ErrorCode=3;
        return;
    end
    
    
    % -> IL CODICE SEGUENTE E' DA CONTROLLARE E VERIFICARE! 
    
    % The Local Machine can be connetted with Remote Computer.
    % Now we verify if user name and password are correct and if remote
    % drive and remote folder exist on the remote computer and it is
    % possible to exchange data with them.
    
    si=[];
    de=[];
    
    [si de]=system(['psinfo \\', DataInput.PcName, ' -u ',DataInput.user, ' -p ',DataInput.passwd ]);
    
    if si<0
        % It is possible to be connected to the remote computer but it is not usable because the user
        % name and/or password is/are incorrect.
        
        ErrorCodeComputer=4;
        return;
    else 
        % Username and Password are correct!         
    end
    
    % Now we verify if it possible to exchange data with the remote
    % computer:
    
    
    fid = fopen('Tracing.txt','w+');
    fclose (fid);

    % ATTENZIONE: verificare perche sembra funzionare anche se il RemoteFolder non
    % esiste.
    
    Status=movefile('Tracing.txt', ['\\',DataInput.PcName,'\',DataInput.RemoteDrive,'$\',DataInput.RemoteFolder]);
    
    if Status==1   
        % Remote Drive/Folder exist on Remote computer and
        % it is possible to exchange data with him.
    else
        
        % Move file error!
        ErrorCodeComputer=5;
        return;
    end
    
    % At this point we can to analyze the remote computer hardware.
    
    
    RealNumCPU=-1;
    RealNumCPU=GiveCPUnumber(de);
    
    % Trasforming the input data provided in a form [n1:n2] in a single numerical
    % value.
    
    
    DataInput.NumCPU=length(DataInput.NumCPU);
    
    if DataInput.NumCPU  == RealNumCPU
        % It is Ok!
    end
    
    if DataInput.NumCPU > RealNumCPU
        ErrorCode=2.1;
        
    end
    if DataInput.NumCPU < RealNumCPU
        ErrorCode=2.2;
    end
    
end   


