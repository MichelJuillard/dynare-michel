function [nCPU]= GiveCPUnumber (ComputerInformations)
% DESCRIPTION
% This function return the CPUs or cores numer avaiable
% on the computer used for run parallel code.
%
% INPUTS
% an array contained several fields that describe the hardaware 
% software enviroments of a generic computer.
%    
% OUTPUTS
% The CPUs or Cores numbers of computer. 
%        
% SPECIAL REQUIREMENTS
%    none

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


nCPU=-1;

OffSet=27;

SringPosition=strfind(ComputerInformations, 'Processors:');
nCPU=ComputerInformations(SringPosition+OffSet);

% We check if there are Processors/Cores more than 9.


t0=ComputerInformations(SringPosition+OffSet+1);
t1=str2num(t0);
t1=isempty(t1);

% if t1 is 0 the machine have more than 9 CPU.

if t1==0
    nCPU=strcat(nCPU,t0);    
end

nCPU=str2num(nCPU);

return
