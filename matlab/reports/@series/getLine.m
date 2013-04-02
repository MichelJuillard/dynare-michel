function o = getLine(o, xrange)
%function o = getLine(o, xrange)
% Create the series
%
% INPUTS
%   o       [series]    series object
%   xrange  [dynDates]  range of x values for line
%
% OUTPUTS
%   o       [series]    series object
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

%% Validate options provided by user
assert(~isempty(o.data) && isa(o.data, 'dynSeries'), ['@series.getLine: must ' ...
                    'provide data as a dynSeries']);

% Line
assert(ischar(o.color), '@series.getLine: color must be a string');
valid_line_style = {'none', '-', '--', ':', '-.'};
assert(any(strcmp(o.line_style, valid_line_style)), ...
       ['@series.getLine: line_style must be one of ' strjoin(valid_line_style, ' ')]);
assert(isfloat(o.line_width), ['@series.getLine: line_width must be a ' ...
                    'positive number']);

% Graph_Marker
valid_graph_marker = {'+', 'o', '*', '.', 'x', 's', 'square', 'd', 'diamond', ...
                '^', 'v', '>', '<', 'p', 'pentagram', 'h', 'hexagram', ...
                'none'};
assert(isempty(o.graph_marker) || any(strcmp(o.graph_marker, valid_graph_marker)), ...
       ['@series.getLine: graph_marker must be one of ' strjoin(valid_graph_marker)]);

assert(ischar(o.graph_marker_edge_color), '@series.getLine: graph_marker_edge_color must be a string');
assert(ischar(o.graph_marker_face_color), '@series.getLine: graph_marker_face_color must be a string');
assert(isfloat(o.graph_marker_size), ['@series.getLine: graph_marker_size must be a ' ...
                    'positive number']);

% Marker & Line
assert(~(strcmp(o.line_style, 'none') && isempty(o.graph_marker)), ['@series.getLine: ' ...
                    'you must provide at least one of line_style and graph_marker']);

% Validate xrange
assert(isempty(xrange) || isa(xrange, 'dynDates'));

%%
if isempty(xrange) || xrange == o.data.time
    ds = o.data;
else
    ds = o.data(xrange);
end

opt = {'XData', 1:length(ds.data)};
opt = {opt{:}, 'YData', ds.data};

opt = {opt{:}, 'Color', o.color};
opt = {opt{:}, 'LineStyle', o.line_style};
opt = {opt{:}, 'LineWidth', o.line_width};

if ~isempty(o.graph_marker)
    opt = {opt{:}, 'Marker', o.graph_marker};
    opt = {opt{:}, 'MarkerSize', o.graph_marker_size};
    opt = {opt{:}, 'MarkerEdgeColor', o.graph_marker_edge_color};
    opt = {opt{:}, 'MarkerFaceColor', o.graph_marker_face_color};
end

line(opt{:});
end
