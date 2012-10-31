% matlab script for testing matlab routines.

% Copyright (C) 2011 Dynare Team
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

list_of_routines = cell(1);

list_of_routines(1,1) = {'../../matlab/@dynSeries'};         % path
list_of_routines(1,2) = {'dynSeries'};                    % name

% list_of_routines = [list_of_routines ; { '../../matlab/kalman/likelihood' , 'kalman_filter' } ];



% Temporarely add ../matlab/utilities/tests to the path.
addpath('../../matlab/utilities/tests')

% Copy output to a file.
dirlog = datestr(now, 'dddd-mmm-yy');
if ~exist(dirlog)
    mkdir(dirlog);
end
diary([dirlog '/errors.log'])

% Run tests.
n = size(list_of_routines,1);
c = zeros(n,1);

for i=1:n
    c(i) = mtest(list_of_routines{1,2},list_of_routines{1,1});
end

if all(c)
    disp('Nothing to report!')
end

diary off
