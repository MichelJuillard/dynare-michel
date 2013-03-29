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
assert(~isempty(o.data) && isa(o.data, 'dynSeries'), ['@series.series: must ' ...
                    'provide data as a dynSeries']);

% Line
assert(ischar(o.color), '@series.series: color must be a string');
valid_line_style = {'none', '-', '--', ':', '-.'};
assert(any(strcmp(o.line_style, valid_line_style)), ...
       ['@series.series: line_style must be one of ' strjoin(valid_line_style, ' ')]);
assert(isfloat(o.line_width), ['@series.series: line_width must be a ' ...
                    'positive number']);

% Marker
valid_marker = {'+', 'o', '*', '.', 'x', 's', 'square', 'd', 'diamond', ...
                '^', 'v', '>', '<', 'p', 'pentagram', 'h', 'hexagram', ...
                'none'};
assert(isempty(o.marker) || any(strcmp(o.marker, valid_marker)), ...
       ['@series.series: marker must be one of ' strjoin(valid_marker)]);

assert(ischar(o.marker_edge_color), '@series.series: marker_edge_color must be a string');
assert(ischar(o.marker_face_color), '@series.series: marker_face_color must be a string');
assert(isfloat(o.marker_size), ['@series.series: marker_size must be a ' ...
                    'positive number']);

% Marker & Line
assert(~(strcmp(o.line_style, 'none') && isempty(o.marker)), ['@series.series: ' ...
                    'you must provide at least one of line_style and marker']);

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

if ~isempty(o.marker)
    opt = {opt{:}, 'Marker', o.marker};
    opt = {opt{:}, 'MarkerSize', o.marker_size};
    opt = {opt{:}, 'MarkerEdgeColor', o.marker_edge_color};
    opt = {opt{:}, 'MarkerFaceColor', o.marker_face_color};
end

line(opt{:});
end
