function [check, info] = mtest(fname, fpath)
% Extract test sections from matlab's routine executes the test and report errors.

% Copyright (C) 2011-2013 Dynare Team
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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

% Default answer (no problem).
check = 1;

% Open the matlab file.
if isempty(fpath)
    % The full path to the matlab routine (with extension) is given.
    fid = fopen(fname,'r');
else
    fid = fopen([fpath '/' fname '.m'],'r');
end

% Read the matlab file.
file = textscan(fid,'%s','delimiter','\n');
file = file{1};

% Close the matlab file.
fclose(fid);

% Locate the test blocks.
b1 = find(strncmp(file,'%@test:',7))+1;
b2 = find(strncmp(file,'%@eof:',6))-1;
nn = length(b2);

if length(b1)-length(b2)
    error('test:: There is a problem with the test blocks definition!')
end

% Initialize the second output if necessary.
if nargout>1
    % First column   name of the tested routine.
    % Second column  number of the unitary test.
    % Third column   status of the unitary test (0 if the test fails, 1 otherwise)
    % Fourth column  details about the failure (vector of integers)
    % Fifth column   elapsed time in seconds (cpu time).
    info = cell(nn,4);
end

% Perform the tests.
for i=1:nn
    if nargout>1
        info(i,1) = {fname};
        info(i,2) = {i};
    end
    % Write the temporary test routine.
    tid = fopen([fname '_test_' int2str(i) '.m'],'w');
    fprintf(tid,['function [T,t,LOG] = ' fname '_test_' int2str(i) '()\n']);
    fprintf(tid,['try\n']);
    for j=b1(i):b2(i)
        str = file{j};
        fprintf(tid,[str(4:end) '\n']);
    end
    fprintf(tid,['LOG = NaN;\n']);
    fprintf(tid,['catch exception\n']);
    fprintf(tid,['LOG = getReport(exception,''extended'');\n']);
    fprintf(tid,['T = NaN;\n']);
    fprintf(tid,['t = NaN;\n']);
    fprintf(tid,['end\n']);
    fclose(tid);
    % Call the temporary test routine.
    init = cputime;
    [TestFlag,TestDetails,LOG] = feval([fname '_test_' int2str(i)]);
    time = cputime-init;
    if isnan(TestFlag)
        fprintf(['\n'])
        fprintf(['Call to ' fname ' test routine n°' int2str(i) ' failed (' datestr(now) ')!\n'])
        fprintf(['\n'])
        disp(LOG)
        check = 0;
        if nargout>1
            info(i,3) = {0};
        end
        continue
    end
    if ~TestFlag
        if nargout>1
            info(i,3) = {0};
            tmp = ones(length(TestDetails),1);
        end
        fprintf(['Test n°' int2str(i) ' for routine ' fname ' failed (' datestr(now) ')!\n']);
        for j=1:length(TestDetails)
            if ~TestDetails(j)
                if nargout>1
                    tmp(j) = 0;
                end
                fprintf(['Output argument n°' int2str(j) ' didn''t give the expected result.\n']);
            end
        end
        info(i,4) = {tmp};
        check = 0;
    else
        if nargout>1
            info(i,3) = {1};
            info(i,4) = {ones(length(TestDetails),1)};
            info(i,5) = {time};
        end
        delete([fname '_test_' int2str(i) '.m'])
    end
end