function o = page(varargin)
%function o = page(varargin)
% Page Class Constructor
%
% INPUTS
%   0 args => empty page
%   1 arg (page class) => copy object
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
o.paper = '';
o.title = {};
o.titleFormat = {};
o.orientation = '';
o.footnote = {};
o.sections = sections();

if nargin == 1
    assert(isa(varargin{1}, 'page'), ['@page.page: with one arg to Page ' ...
                        'constructor, you must pass a page object']);
    o = varargin{1};
    return;
elseif nargin > 1
    if round(nargin/2) ~= nargin/2
        error(['@page.page: options must be supplied in name/value ' ...
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
            error('@page.page: %s is not a recognized option.', pair{1});
        end
    end
end

% Check options provided by user
if ischar(o.title)
    o.title = {o.title};
end
if ischar(o.titleFormat)
    o.titleFormat = {o.titleFormat};
end
assert(iscellstr(o.title), ...
       '@page.page: title must be a cell array of strings');
assert(iscellstr(o.titleFormat), ...
       '@page.page: titleFormat must be a cell array of strings');
assert(length(o.title)==length(o.titleFormat), ...
       '@page.page: title and titleFormat must be of the same length');

valid_paper = {'a4', 'letter'};
assert(any(strcmp(o.paper, valid_paper)), ...
       ['@page.page: paper must be one of ' strjoin(valid_paper, ' ')]);

valid_orientation = {'portrait', 'landscape'};
assert(any(strcmp(o.orientation, valid_orientation)), ...
       ['@page.page: orientation must be one of ' strjoin(valid_orientation, ' ')]);

if ischar(o.footnote)
    o.footnote = {o.footnote};
end
assert(iscellstr(o.footnote), ...
       '@page.page: footnote must be a cell array of string(s)');

% Create page object
o = class(o, 'page');
end
