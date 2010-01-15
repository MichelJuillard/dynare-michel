function fMessageStatus(prtfrc, njob, waitbarString, waitbarTitle, Parallel, MasterName, DyMo)
% In a parallelization context, this function is launched on slave
% machines, and acts as a message passing device for the master machine.
% prtfrc [double], fraction of iteration performed
% njob [int], index number of this CPU among all CPUs in the
%                       cluster
% waitbarString [char], running message string to be displayed in the monitor window on master machine 
% waitbarTitle [char], title to be displayed in the monitor window on master machine
% Parallel [struct], options_.parallel(ThisMatlab), i.e. the paralle settings for this slave machine in the cluster 
% MasterName [char], IP address or PC name of the master 
% DyMo [char], working directory of the master machine

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

global funcName

if nargin<5,
    Parallel.Local=1;
end

save(['comp_status_',funcName,int2str(njob),'.mat'],'prtfrc','njob','waitbarString','waitbarTitle');
if Parallel.Local==0,
    if isunix || (~matlab_ver_less_than('7.4') && ismac),
        system(['scp comp_status_',funcName,int2str(njob),'.mat ',Parallel.user,'@',MasterName,':',DyMo]);
    else
        copyfile(['comp_status_',funcName,int2str(njob),'.mat'],['\\',MasterName,'\',DyMo(1),'$\',DyMo(4:end),'\']);
    end
end
