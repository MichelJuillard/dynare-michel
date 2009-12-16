function fMessageStatus(prtfrc, njob, waitbarString, waitbarTitle, Parallel, MasterName, DyMo)

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
    if isunix,
        system(['scp comp_status_',funcName,int2str(njob),'.mat ',Parallel.user,'@',MasterName,':',DyMo]);
    else
        copyfile(['comp_status_',funcName,int2str(njob),'.mat'],['\\',MasterName,'\',DyMo(1),'$\',DyMo(4:end),'\']);
    end
end
