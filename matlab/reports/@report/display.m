function display(r)
%function display(r)
% Display a Report object
%
% INPUTS
%   none
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

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

disp(' ');
disp([inputname(1) '.title = ']);
disp(' ');
disp(['     ''' r.title '''']);
disp(' ')
disp([inputname(1) '.orientation = ']);
disp(' ');
disp(['     ''' r.orientation '''']);
disp(' ')
disp([inputname(1) '.numPages() = ']);
disp(' ');
disp(['     ' num2str(numPages(r))]);
disp(' ');
disp([inputname(1) '.pages = ']);
disp(' ');
disp(r.pages.getPages());
end