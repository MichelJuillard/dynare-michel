function closeSlave(Parallel),
% In a parallelc context, this utility closes all remote matlab instances
% called by masteParallelMan (which leaves open remote matlab instances)

% Copyright (C) 2010 Dynare Team
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

delete( 'slaveParallel_input*.mat');
for indPC=1:length(Parallel),
    if Parallel(indPC).Local==0,
        if isunix || (~matlab_ver_less_than('7.4') && ismac),
            system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' rm -fr ',Parallel(indPC).RemoteFolder,'/slaveParallel_input*.mat']);
        else
            mydelete('slaveParallel_input*.mat',['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\']);
        end
    end
end
