function [nCPU, totCPU, nBlockPerCPU] = distributeJobs(Parallel, fBlock, nBlock)

% Determine the total number of available CPUs, and the number of threads to run on each CPU

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

totCPU=0;
for j=1:length(Parallel),
    nCPU(j)=length(Parallel(j).NumCPU);
    totCPU=totCPU+nCPU(j);
end

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
