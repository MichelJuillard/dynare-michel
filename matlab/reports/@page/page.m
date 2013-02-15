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
o.caption = '';
o.orientation = 'portrait';
o.sections = sections();

if nargin == 1
    assert(isa(varargin{1}, 'page'),['With one arg to Page constructor, ' ...
                        'you must pass a page object']);
    o = varargin{1};
    return;
elseif nargin > 1
    if round(nargin/2) ~= nargin/2
        error(['Options to Page constructor must be supplied in name/value ' ...
               'pairs.']);
    end

    optNames = lower(fieldnames(o));

    % overwrite default values
    for pair = reshape(varargin, 2, [])
        field = lower(pair{1});
        if any(strmatch(field, optNames, 'exact'))
            o.(field) = pair{2};
        else
            error('%s is not a recognized option to the Page constructor.', ...
                  field);
        end
    end
end

% Create page object
o = class(o, 'page');
end
