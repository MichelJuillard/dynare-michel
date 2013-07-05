function report = run_unitary_tests_in_directory(dirname,savereport,printreport,sendreport)

% Copyright (C) 2013 Dynare Team
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

system('git show --pretty=format:"Last commit %H by %an, %ar %n-> %s" HEAD > git.info');
system('git rev-parse HEAD > git.last-commit-hash');

fid = fopen('git.info');
gitinfo = fgetl(fid);
gitinfo = char(gitinfo,fgetl(fid));
fclose(fid);

fid = fopen('git.last-commit-hash');
gitlastcommithash = fgetl(fid);
fclose(fid);

matlabverion = version;
platform = computer;

listoffiles = get_directory_description(dirname);

diary(['report-' gitlastcommithash '.log'] )
[report, time] = run_unitary_tests(listoffiles)
diary off

if nargin>1 && savereport>0
    save(['report-' gitlastcommithash '.mat'],'report','time','gitinfo','gitlastcommithash','matlabverion','platform');
end

if nargin>2
    build_report_summary(['report-' gitlastcommithash '.mat'], printreport, sendreport);
end