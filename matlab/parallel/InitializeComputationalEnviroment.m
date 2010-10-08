function [ErrorCode] = InitializeComputationalEnviroment(DataInput)

% PARALLEL CONTEXT
% In a parallel context, this function is used to Initialize the computational enviroment according with
% the user request.
% If no error happen the function return 0.
%
% INPUT/OUTPUT description:
%
%
%
%
% The variable ErrorCode is initialized at 0. If there are non problems 
% the ErrorCode is unchanged, in the others cases is set equak to 1, 2 , ... The values
% table is below.
%
%
%   Table for ErrorCode Values.
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

% Invoke masterParallel with 8 arguments and the last equal to 1. With this shape
% for input data, masterParallel only create a new directory for remote
% computation. The name of this directory is time depending. For local
% parallel computations with Strategy == 1 delete the traces (if exists) of
% previous computations. 
delete(['P_slave_*End.txt'])
masterParallel(DataInput.parallel,[],[],[],[],[],[],DataInput.parallel_info,1);

return