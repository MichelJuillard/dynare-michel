function dynareParallelMkDir(PRCDir,Parallel)
% PARALLEL CONTEXT
% In a parallel context, this is a specialized version of rmdir() function.
%
%
% INPUT/OUTPUT description:
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



if nargin ==0,
    disp('dynareParallelMkDir(dirname,Parallel)')
    return
end

for indPC=1:length(Parallel)
    if Parallel(indPC).Local==0,
        if isunix || (~matlab_ver_less_than('7.4') && ismac),
            system(['ssh ',Parallel(indPC).user,'@',Parallel(indPC).PcName,' mkdir -p ',Parallel(indPC).RemoteFolder,'/',PRCDir])
        else
            [NonServeS NonServeD]=mkdir(['\\',Parallel(indPC).PcName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteFolder,'\',PRCDir]);
        end
    end
end

return