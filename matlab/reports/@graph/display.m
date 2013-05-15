function display(o)
%function display(o)
% Display a Graph object
%
% INPUTS
%   o   [graph] graph object
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

name = 'report.page.section.graph';
disp(' ');
disp([name '.title = ']);
disp(' ');
disp(['     ''' o.title '''']);

disp(' ');
disp([name '.xlabel = ']);
disp(' ');
disp(['     ''' o.xlabel '''']);

disp(' ');
disp([name '.ylabel = ']);
disp(' ');
disp(['     ''' o.ylabel '''']);

disp(' ');
disp([name '.footnote = ']);
disp(' ');
disp(['     ''' o.footnote '''']);

disp(' ');
disp([name '.figname = ']);
disp(' ');
disp(['     ''' o.figname '''']);

disp(' ');
disp([name '.legend = ']);
disp(' ');
disp(o.legend);

disp(' ');
disp([name '.shade = ']);
disp(' ');
disp(o.shade);

disp(' ');
disp([name '.seriesElements = ']);
disp(' ');
o.seriesElements.getSeriesElements()
end