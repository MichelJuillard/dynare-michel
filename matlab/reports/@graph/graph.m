function o = graph(varargin)
%function o = graph(varargin)
% Graph Class Constructor
%
% INPUTS
%   varargin        0 args  : empty graph object
%                   1 arg   : must be graph object (return a copy of arg)
%                   > 1 args: option/value pairs (see structure below for
%                   options)
%
% OUTPUTS
%   o   [graph] graph object
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

o = struct;

o.config = '';

o.title = '';
o.ylabel = '';
o.xlabel = '';
o.footnote = '';

o.figname = '';
o.data = '';
o.seriestouse = '';
o.shade = '';
o.xrange = '';
o.yrange = '';

o.grid = true;

o.legend = false;
o.legend_location = 'SouthEast';
o.legend_orientation = 'horizontal';
o.legend_font_size = 8;

if nargin == 1
    assert(isa(varargin{1}, 'graph'),['@graph.graph: with one arg you ' ...
                        'must pass a graph object']);
    o = varargin{1};
    return;
elseif nargin > 1
    if round(nargin/2) ~= nargin/2
        error(['@graph.graph: options must be supplied in name/value ' ...
               'pairs.']);
    end

    optNames = lower(fieldnames(o));

    % overwrite default values
    for pair = reshape(varargin, 2, [])
        field = lower(pair{1});
        if any(strmatch(field, optNames, 'exact'))
            o.(field) = pair{2};
        else
            error('@graph.graph: %s is not a recognized option.', field);
        end
    end
end

% Check options provided by user
assert(ischar(o.title), '@graph.graph: title must be a string');
assert(ischar(o.footnote), '@graph.graph: footnote must be a string');
assert(ischar(o.config), '@graph.graph: config file must be a string');
assert(ischar(o.xlabel), '@graph.graph: xlabel file must be a string');
assert(ischar(o.ylabel), '@graph.graph: ylabel file must be a string');
assert(ischar(o.figname), '@graph.graph: figname must be a string');
assert(islogical(o.grid), '@graph.graph: grid must be either true or false');
assert(islogical(o.legend), '@graph.graph: legend must be either true or false');
assert(isint(o.legend_font_size), '@graph.graph: legend_font_size must be an integer');
valid_legend_locations = ...
    {'North', 'South', 'East', 'West', ...
     'NorthEast', 'SouthEast', 'NorthWest', 'SouthWest', ...
     'NorthOutside', 'SouthOutside', 'EastOutside', 'WestOutside', ...
     'NorthEastOutside', 'SouthEastOutside', 'NorthWestOutside', 'SouthWestOutside', ...
     'Best', 'BestOutside', ...
    };
assert(any(strcmp(o.legend_location, valid_legend_locations)), ...
       ['@graph.graph: legend_location must be one of ' strjoin(valid_legend_locations, ' ')]);

valid_legend_orientations = {'vertical', 'horizontal'};
assert(any(strcmp(o.legend_orientation, valid_legend_orientations)), ...
       ['@graph.graph: legend_orientation must be one of ' strjoin(valid_legend_orientations, ' ')]);

assert(isempty(o.shade) || (iscell(o.shade) && length(o.shade) == 2 && ...
                            ischar(o.shade{1}) && ischar(o.shade{2})), ...
       ['@graph.graph: yrange is specified as ''{''1999q1'',''1999q2''}''.']);
assert(isempty(o.xrange) || (iscell(o.xrange) && length(o.xrange) == 2 && ...
                             ischar(o.xrange{1}) && ischar(o.xrange{2})), ...
       ['@graph.graph: xrange is specified as ''{''1999q1'',''1999q2''}''.']);
assert(isempty(o.yrange) || (iscell(o.yrange) && length(o.yrange) == 2 && ...
                             ischar(o.yrange{1}) && ischar(o.yrange{2})), ...
       ['@graph.graph: yrange is specified as ''{''1999q1'',''1999q2''}''.']);

assert(~isempty(o.data), '@graph.graph: must provide data');
msg = ['@graph.graph: data must either be a dynSeries or a cell array of ' ...
       'dynSeries'];
if length(o.data) == 1
    assert(isa(o.data, 'dynSeries'), msg);
else
    assert(iscell(o.data), msg);
    for i=1:length(o.data)
        assert(isa(o.data{i}, 'dynSeries'), msg);
    end
end

msg = ['@graph.graph: series to use must be a cell array of string(s)'];
if ~isempty(o.seriestouse)
    assert(iscell(o.seriestouse), msg);
    for i=1:length(o.seriestouse)
        assert(ischar(o.seriestouse{i}), msg);
    end
end

% Create graph object
o = class(o, 'graph');
end