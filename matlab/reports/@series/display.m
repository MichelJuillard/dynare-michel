function display(o)
%function display(o)
% Display a Series object
%
% INPUTS
%   o   [series] series object
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

name = 'report.page.section.graph.series';
disp(' ');
disp([name '.data = ']);
disp(' ');
display(o.data);

disp(' ');
disp([name '.color = ']);
disp(' ');
disp(['     ''' o.color '''']);

disp(' ');
disp([name '.line_style = ']);
disp(' ');
disp(['     ''' o.line_style '''']);

disp(' ');
disp([name '.line_width = ']);
disp(' ');
disp(['     ''' o.line_width '''']);

disp(' ');
disp([name '.marker = ']);
disp(' ');
disp(['     ''' o.marker '''']);

disp(' ');
disp([name '.marker_edge_color = ']);
disp(' ');
disp(['     ''' o.marker_edge_color '''']);

disp(' ');
disp([name '.marker_face_color = ']);
disp(' ');
disp(['     ''' o.marker_face_color '''']);

disp(' ');
disp([name '.marker_size = ']);
disp(' ');
disp(['     ''' o.marker_size '''']);
end