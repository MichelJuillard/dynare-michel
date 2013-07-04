function report = run_unitary_tests_in_directory(dirname,savereport)

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

listoffiles = get_directory_description(dirname);
report = run_unitary_tests(listoffiles);

if nargin>1 && savereport>0
    d = clock;
    save(['report-' num2str(d(1)) '-' num2str(d(2)) '-' num2str(d(3)) '-' num2str(d(4)) '-' num2str(d(5)) '-'  num2str(round(d(6))) '.mat' ],'report');
end