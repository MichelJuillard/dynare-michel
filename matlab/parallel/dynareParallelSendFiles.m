function dynareParallelSendFiles(NamFileInput,PRCDir,Parallel)
% PARALLEL CONTEXT
% In a parallel context, this is a specialized mono-directional (Local to Remote) version of copy()
% function.
%
%
% INPUTS
%  o NamFileInput   []   ...
%  o PRCDir         []   ... 
%  o Parallel       []   ...  
%
%  OUTPUTS
%  None
%
%
%
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


if ischar(NamFileInput),
    for j=1:size(NamFileInput,1),
        NamFile(j,:)={'',deblank(NamFileInput(j,:))};
    end
    NamFileInput = NamFile;
end

for indPC=1:length(Parallel),
    if Parallel(indPC).Local==0,
        if isunix || (~matlab_ver_less_than('7.4') && ismac),
            for jfil=1:size(NamFileInput,1),
                if ~isempty(NamFileInput{jfil,1})
                    system(['ssh ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,' mkdir -p ',Parallel(indPC).RemoteDirectory,'/',PRCDir,'/',NamFileInput{jfil,1}])
                end
                system(['scp ',NamFileInput{jfil,1},NamFileInput{jfil,2},' ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,':',Parallel(indPC).RemoteDirectory,'/',PRCDir,'/',NamFileInput{jfil,1}]);
            end
        else
            for jfil=1:size(NamFileInput,1)
                if ~isempty(NamFileInput{jfil,1})
                    if isempty(dynareParallelDir(NamFileInput{jfil,1},PRCDir,Parallel(indPC)))
                        mkdir(['\\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\',NamFileInput{jfil,1}]);
                    end
                end
                copyfile([NamFileInput{jfil,1},NamFileInput{jfil,2}],['\\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\',NamFileInput{jfil,1}])
            end
        end
    end
end
