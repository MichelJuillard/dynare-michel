function [ErrorCode] = AnalyseComputationalEnviroment(DataInput)


% PARALLEL CONTEXT
% In a parallel context, this function is used to check the user computational request.
% If no error happen the function return 0.
%
% INPUT/OUTPUT description:
%
% DataInput is the strcture option_.parallel, with the follow fields:
%
%           Local         Define the computation place: 1 is on local machine, 0 remote
%     ComputerName        Intuitive: contain the computer name.
%           CPUnbr        Intuitive: contain the CPU number.
%         UserName        Intuitive: contain the use name for the ComputerName.
%           passwd        Intuitive: contain the password for the user name in ComputerName.
%      RemoteDrive        Drive used for Local/Remote computation (data exchange, etc) must be contain 'RemoteFolder'.
%     RemoteDirectory     Folder in RemoteDrive used for Local/Remote computation.
%     MatlabOctavePath    []
%         DynarePath      []
%
%   This information is typed by the user using the *.mod file, 
%   the goal of this function is to check if it correct.
%
%
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
%               have s>d!
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
%   Value 7:    The valu user and/or passwd are incorret or the user have
%               no permission to execute a Matlab section.
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


% The function is composed by two main blocks, determined by the 'Local'
% variable. 

if ((DataInput.Local == 0) |(DataInput.Local == 1))
    % Continue it is Ok!
else
    ErrorCode=1;
    return

end


%%%%%%%%%%  Local Machine   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this case we need to check only the variable 'CPUnbr'. 

% We run the parallel code on local computer, so the others fields are automatically 
% fixed by Dynare. Then the user can also fill them with wrong values.

if (DataInput.Local == 1)
    
    yn=isempty(DataInput.CPUnbr);
    
    if yn==1
        ErrorCode=2;
        return
    end
    
    % We look for the information on local computer hardware.
    
    si0=[];
    de0=[];
    
    [si0 de0]=system(['psinfo \\']);
    
    RealCPUnbr=-1;
    RealCPUnbr=GiveCPUnumber(de0);
    
    % Trasforming the input data provided in a form [n1:n2] in a single numerical
    % value.   
    
    DataInput.CPUnbr=length(DataInput.CPUnbr);
    
    if DataInput.CPUnbr  == RealCPUnbr
        % It is Ok!
    end
    
    if DataInput.CPUnbr > RealCPUnbr
        ErrorCode=2.1;
        
    end
    if DataInput.CPUnbr < RealCPUnbr
        ErrorCode=2.2;
    end    
end   

%%%%%%%%%%  Remote Machine   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this case we need more sophisticated check. 

 
if (DataInput.Local == 0)
    
   si1=[];
   de1=[]; 
    
  [si1 de1]=system(['ping ', DataInput.ComputerName]);
    
    if si1==1 
        % It is impossiblie to be connected to the
        % remote computer.
        
        ErrorCode=3;
        return;
    end
    
    % The Local Machine can be connetted with Remote Computer.
    
    % Now we verify if user name and password are correct and if remote
    % drive and remote folder exist on the remote computer and it is
    % possible to exchange data with them.
   
     
    if (isempty(DataInput.UserName)) || (isempty(DataInput.Password))
      
      % The fields user name and/or password are/is empty!
      
      ErrorCode=4;
      return
    
    end
    
   
    % Now we very if RemoteDrive and/or RemoteDirectory exist on remote
    % computer
    
    StartPwd=pwd;
    
    if (isempty(DataInput.RemoteDirectory))
        disp 'Remote Drive and/or Remote Folder not exist!';
        
        ErrorCode=5;
        return
    end
    
    
    try
        cd(['\\',DataInput.ComputerName,'\',DataInput.RemoteDrive,'$\',DataInput.RemoteDirectory]);
    catch
       
        cd ([StartPwd]);
        disp 'Remote Drive and/or Remote Folder not exist!';
        
        ErrorCode=5;
        return
    
    end

 
    
   cd ([StartPwd]);
    
 
   % Now we verify if it possible to exchange data with the remote
   % computer:
    
    fT = fopen('WriteTest.txt', 'w+');
    fclose(fT);
   
    Status=copyfile('WriteTest.txt', ['\\',DataInput.ComputerName,'\',DataInput.RemoteDrive,'$\',DataInput.RemoteDirectory]);
    
    
    if Status==1   
        
        % Remote Drive/Folder exist on Remote computer and
        % it is possible to exchange data with him.
        
        %Delete Traces (local and remotely)
        delete ('WriteTest.txt')
        delete(['\\',DataInput.ComputerName,'\',DataInput.RemoteDrive,'$\',DataInput.RemoteDirectory,'\WriteTest.txt']);
    else
        
        ErrorCode=6;
        return;
    end
   
    
   % Now we verify if it is possible execute a matlab section on remote
   % machine when the user is .UserName with password .Password
   
    si2=[];
    de2=[];
   
    [si2 de2]=system(['start /B /WAIT psexec \\',DataInput.ComputerName,' -e -u ',DataInput.UserName,' -p ',DataInput.Password,' -W ',DataInput.RemoteDrive,':\',DataInput.RemoteDirectory, ' -low  matlab -nosplash -nodesktop -minimize -r Tracing']);
     
    NoError='error code 0';
    
    StrError= findstr(NoError,de2);
    
   
    if isempty (StrError)
      
       % Bad user and/or passwd!
        ErrorCode=7;
        return;
       
    else 
        
        % No error it is possible execute a matlab section on remote
        % machine when the user is .UserName with password .Password
    end
    
   
    % At this point we can to analyze the remote computer hardware.
    
    RealCPUnbr=-1;
    
    
    [si0 de0]=system(['psinfo \\']);
    RealCPUnbr=GiveCPUnumber(de0);
    
    
    % Trasforming the input data provided in a form [n1:n2] in a single numerical
    % value.
    
    
    DataInput.CPUnbr=length(DataInput.CPUnbr);
    
    if DataInput.CPUnbr  == RealCPUnbr
        % It is Ok!
    end
    
    if DataInput.CPUnbr > RealCPUnbr
        ErrorCode=2.1;
        
    end
    if DataInput.CPUnbr < RealCPUnbr
        ErrorCode=2.2;
    end
    
end   


