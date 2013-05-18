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

o.seriesElements = seriesElements();

o.title = '';
o.ylabel = '';
o.xlabel = '';

o.figname = '';
o.data = '';
o.seriesToUse = '';
o.xrange = '';
o.yrange = '';

o.shade = '';
o.shadeColor = 'green';
o.shadeOpacity = .2;

o.showGrid = true;

o.showLegend = false;
o.showLegendBox = false;
o.legendLocation = 'SouthEast';
o.legendOrientation = 'horizontal';
o.legendFontSize = 8;

o.showZeroline = false;

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

    optNames = fieldnames(o);

    % overwrite default values
    for pair = reshape(varargin, 2, [])
        ind = strmatch(lower(pair{1}), lower(optNames), 'exact');
        assert(isempty(ind) || length(ind) == 1);
        if ~isempty(ind)
            o.(optNames{ind}) = pair{2};
        else
            error('@graph.graph: %s is not a recognized option.', pair{1});
        end
    end
end

% Check options provided by user
assert(ischar(o.title), '@graph.graph: title must be a string');
assert(ischar(o.xlabel), '@graph.graph: xlabel file must be a string');
assert(ischar(o.ylabel), '@graph.graph: ylabel file must be a string');
assert(ischar(o.figname), '@graph.graph: figname must be a string');
assert(islogical(o.showGrid), '@graph.graph: showGrid must be either true or false');
assert(islogical(o.showLegend), '@graph.graph: showLegend must be either true or false');
assert(islogical(o.showLegendBox), '@graph.graph: showLegendBox must be either true or false');
assert(isint(o.legendFontSize), '@graph.graph: legendFontSize must be an integer');
assert(islogical(o.showZeroline), '@graph.graph: showZeroline must be either true or false');
assert(ischar(o.shadeColor), '@graph.graph: shadeColor must be a string');
assert(isfloat(o.shadeOpacity) && length(o.shadeOpacity)==1 && ...
       o.shadeOpacity >= 0 && o.shadeOpacity <= 1, ...
       '@graph.graph: o.shadeOpacity must be a real in [0 1]');
valid_legend_locations = ...
    {'North', 'South', 'East', 'West', ...
     'NorthEast', 'SouthEast', 'NorthWest', 'SouthWest', ...
     'NorthOutside', 'SouthOutside', 'EastOutside', 'WestOutside', ...
     'NorthEastOutside', 'SouthEastOutside', 'NorthWestOutside', 'SouthWestOutside', ...
     'Best', 'BestOutside', ...
    };
assert(any(strcmp(o.legendLocation, valid_legend_locations)), ...
       ['@graph.graph: legendLocation must be one of ' strjoin(valid_legend_locations, ' ')]);

valid_legend_orientations = {'vertical', 'horizontal'};
assert(any(strcmp(o.legendOrientation, valid_legend_orientations)), ...
       ['@graph.graph: legendOrientation must be one of ' strjoin(valid_legend_orientations, ' ')]);

assert(isempty(o.shade) || (isa(o.shade, 'dynDates') && o.shade.ndat >= 2), ...
       ['@graph.graph: shade is specified as a dynDates range, e.g. ' ...
        '''dynDates(''1999q1''):dynDates(''1999q3'')''.']);
assert(isempty(o.xrange) || (isa(o.xrange, 'dynDates') && o.xrange.ndat >= 2), ...
       ['@graph.graph: xrange is specified as a dynDates range, e.g. ' ...
        '''dynDates(''1999q1''):dynDates(''1999q3'')''.']);
assert(isempty(o.yrange) || (isfloat(o.yrange) && length(o.yrange) == 2 && ...
                             o.yrange(1) < o.yrange(2)), ...
       ['@graph.graph: yrange is specified an array with two float entries, ' ...
        'the lower bound and upper bound.']);
assert(isempty(o.data) || isa(o.data, 'dynSeries'), ['@graph.graph: data must ' ...
                    'be a dynSeries']);
assert(isempty(o.seriesToUse) || iscellstr(o.seriesToUse), ['@graph.graph: ' ...
                    'series to use must be a cell array of string(s)']);

% using o.seriesToUse, create series objects and put them in o.seriesElements
if ~isempty(o.data)
    if isempty(o.seriesToUse)
        for i=1:o.data.vobs
            o.seriesElements = o.seriesElements.addSeries('data', o.data{o.data.name{i}});
        end
    else
        for i=1:length(o.seriesToUse)
            o.seriesElements = o.seriesElements.addSeries('data', o.data{o.seriesToUse{i}});
        end
    end
end
o = rmfield(o, 'seriesToUse');
o = rmfield(o, 'data');

% Create graph object
o = class(o, 'graph');
end