function closeSlave(Parallel,TmpFolder),
% PARALLEL CONTEXT
% In parallel context, this utility closes all remote matlab instances
% called by masterParallel with strategy (1) i.e. always open (which leaves
% open remote matlab instances).
%
% INPUTS
%  o Parallel [struct vector]   copy of options_.parallel.
%
% OUTPUTS
%   None
%
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

for indPC=1:length(Parallel),
    if (Parallel(indPC).Local==0),
        dynareParallelDelete( 'slaveParallel_input*.mat',TmpFolder,Parallel(indPC));
    else
        delete( 'slaveParallel_input*.mat');
        pause(1)
        delete(['slaveParallel_*.log']);
        
    end
end

while(1)
    if isempty(dynareParallelDir(['P_slave_',int2str(j),'End.txt'],TmpFolder,Parallel));
        for indPC=1:length(Parallel),
            if (Parallel(indPC).Local==0),
                dynareParallelRmDir(TmpFolder,Parallel(indPC)),
                
            end
        end
        break
        
    end
end

        