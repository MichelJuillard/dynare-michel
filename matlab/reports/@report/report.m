function r = report(varargin)
%function r = report(varargin)
% Report Class Constructor
%
% INPUTS
%   0 args => no title, portrait orientation
%   1 arg (report class) => copy class
%   1 arg (not report class) => title
%   2 args (1st not report class) => title, orientation
%   3 args (1st not report class) => title, orientation, configuraiton
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

% default values
r = struct;
r.title = '';
r.orientation = 'portrait';
r.pages = pages();
r.config = '';

% Check arguments
if nargin == 1
    assert(isa(varargin{1}, 'report') || ischar(varargin{1}), ...
        ['With one arg to report constructor, you must pass either '...
        'a report object or a char.']);
end

if nargin > 1
    assert(~isa(varargin{1}, 'report'), ...
        ['With multiple args to report constructor, first argument ' ...
        'cannot be a report object.']);
    for i=1:nargin
        assert(ischar(varargin{i}), ...
            ['With muliple args to report constructor, all '...
            'arguments must be char.']);
    end
end

% Initialize fields
switch nargin
    case 0
    case 1
        if isa(varargin{1}, 'report')
            r = varargin{1};
            return
        else
            r.title = varargin{1};
        end
    case 2
        r.title = varargin{1};
        r.orientation = validateOrientation(varargin{2});
    case 3
        r.title = varargin{1};
        r.orientation = validateOrientation(varargin{2});
        r.config = varargin{3};
    otherwise
        error('Report constructor: invalid number of arguments');
end

% Create report object
r = class(r, 'report');
end