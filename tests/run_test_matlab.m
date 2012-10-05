% Copyright (C) 2011-2012 Dynare Team
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

top_test_dir = getenv('TOP_TEST_DIR');
addpath(top_test_dir);
addpath([top_test_dir filesep '..' filesep 'matlab']);

% Test Dynare Version
if ~strcmp(dynare_version(), getenv('DYNARE_VERSION'))
  error('Incorrect version of Dynare is being tested')
end

% Test MOD files listed in Makefile.am
[modfile, name] = strtok(getenv('FILESTEM'));
[directory, testfile, ext] = fileparts([top_test_dir '/' modfile]);
cd(directory);

disp('');
disp(['***  TESTING: ' modfile ' ***']);
try
  dynare([testfile ext], 'console')
  testFailed = false;
catch exception
  printMakeCheckMatlabErrMsg(strtok(getenv('FILESTEM')), exception);
  testFailed = true;
end

cd(getenv('TOP_TEST_DIR'));
name = strtok(getenv('FILESTEM'));
fid = fopen([name '.m.trs'], 'w');
if fid < 0
  wd = pwd
  filestep = getenv('FILESTEM')
  error(['ERROR: problem opening file ' name '.m.trs for writing....']);
end
if testFailed
  fprintf(fid,':test-result: FAIL\n');
  fprintf(fid,':number-tests: 1\n');
  fprintf(fid,':number-failed-tests: 1\n');
  fprintf(fid,':list-of-failed-tests: %s\n', [name '.mod']);
  fprintf(fid,':copy-in-global-log: yes\n');
  fprintf(fid,':recheck: yes\n');
  fprintf(fid,':test-global-result: FAIL\n');
else
  fprintf(fid,':test-result: PASS\n');
  fprintf(fid,':number-tests: 1\n');
  fprintf(fid,':number-failed-tests: 0\n');
  fprintf(fid,':list-of-failed-tests: \n');
  fprintf(fid,':copy-in-global-log: no\n');
  fprintf(fid,':recheck: no\n');
  fprintf(fid,':test-global-result: PASS\n');
end
fclose(fid);
exit;
