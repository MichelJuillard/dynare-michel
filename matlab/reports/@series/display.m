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
disp([name '.graphLineColor = ']);
disp(' ');
disp(['     ''' o.graphLineColor '''']);

disp(' ');
disp([name '.graphLineStyle = ']);
disp(' ');
disp(['     ''' o.graphLineStyle '''']);

disp(' ');
disp([name '.graphLineWidth = ']);
disp(' ');
disp(['     ''' o.graphLineWidth '''']);

disp(' ');
disp([name '.graphMarker = ']);
disp(' ');
disp(['     ''' o.graphMarker '''']);

disp(' ');
disp([name '.graphMarkerEdgeColor = ']);
disp(' ');
disp(['     ''' o.graphMarkerEdgeColor '''']);

disp(' ');
disp([name '.graphMarkerFaceColor = ']);
disp(' ');
disp(['     ''' o.graphMarkerFaceColor '''']);

disp(' ');
disp([name '.graphMarkerSize = ']);
disp(' ');
disp(['     ''' o.graphMarkerSize '''']);

disp(' ');
disp([name '.tableAlignRight = ']);
disp(' ');
disp(['     ''' o.tableAlignRight '''']);

disp(' ');
disp([name '.tableNegColor = ']);
disp(' ');
disp(['     ''' o.tableNegColor '''']);

disp(' ');
disp([name '.tablePosColor = ']);
disp(' ');
disp(['     ''' o.tablePosColor '''']);

disp(' ');
disp([name '.tableShowMarkers = ']);
disp(' ');
disp(['     ''' o.tableShowMarkers '''']);

end