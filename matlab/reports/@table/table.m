function o = table(varargin)
%function o = table(varargin)
% Table Class Constructor
%
% INPUTS
%   0 args => empty table
%   1 arg (table class) => copy object
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

o = struct;

o.seriesElements = seriesElements();

o.title = '';
o.titleSize = 'large';

o.showHlines = false;
o.showVlines = false;
o.vlineAfter = '';

o.data = '';
o.seriesToUse = '';
o.range = {};
o.precision = 1;

if nargin == 1
    assert(isa(varargin{1}, 'table'),['With one arg to Table constructor, ' ...
                        'you must pass a table object']);
    o = varargin{1};
    return;
elseif nargin > 1
    if round(nargin/2) ~= nargin/2
        error(['Options to Table constructor must be supplied in name/value ' ...
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
            error('%s is not a recognized option to the Table constructor.', pair{1});
        end
    end
end

% Check options provided by user
assert(ischar(o.title), '@table.table: title must be a string');
assert(islogical(o.showHlines), '@table.table: showHlines must be true or false');
assert(islogical(o.showVlines), '@table.table: showVlines must be true or false');
assert(isint(o.precision), '@table.table: precision must be an int');
assert(isempty(o.range) || (isa(o.range, 'dynDates') && o.range.ndat >= 2), ...
       ['@table.table: range is specified as a dynDates range, e.g. ' ...
        '''dynDates(''1999q1''):dynDates(''1999q3'')''.']);
assert(isempty(o.data) || isa(o.data, 'dynSeries'), ...
       '@table.table: data must be a dynSeries');
assert(isempty(o.seriesToUse) || iscellstr(o.seriesToUse), ...
       '@table.table: seriesToUse must be a cell array of string(s)');
assert(isempty(o.vlineAfter) || isa(o.vlineAfter, 'dynDate'), ...
       '@table.table: vlineAfter must be a dynDate');
if o.showVlines
    o.vlineAfter = '';
end
valid_title_sizes = {'Huge', 'huge', 'LARGE', 'Large', 'large', 'normalsize', ...
                    'small', 'footnotesize', 'scriptsize', 'tiny'};
assert(any(strcmp(o.titleSize, valid_title_sizes)), ...
       ['@table.table: titleSize must be one of ' strjoin(valid_title_sizes, ' ')]);

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

% Create table object
o = class(o, 'table');
end