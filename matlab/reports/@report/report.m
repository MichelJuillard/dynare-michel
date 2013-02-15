function o = report(varargin)
%function o = report(varargin)
% Report Class Constructor
%
% INPUTS
%   1 report class object => make a copy
%   Otherwise, option/value pairs (see structure below for options)
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
o = struct;
o.title = '';
o.orientation = 'portrait';
o.paper = 'a4';
o.pages = pages();
o.filename = 'report.tex';
o.config = '';

if nargin == 1
    assert(isa(varargin{1}, 'report'),['With one arg to Report constructor, ' ...
                        'you must pass a report object']);
    r = varargin{1};
    return;
elseif nargin > 1
    if round(nargin/2) ~= nargin/2
        error(['Options to Report constructor must be supplied in name/value ' ...
               'pairs.']);
    end

    optNames = lower(fieldnames(o));

    % overwrite default values
    for pair = reshape(varargin, 2, [])
        field = lower(pair{1});
        if any(strmatch(field, optNames, 'exact'))
            if strcmp(field, 'orientation')
                validateOrientation(pair{2});
            elseif strcmp(field, 'paper')
                validatePaper(pair{2});
            end
            o.(field) = pair{2};
        else
            error('%s is not a recognized option to the Report constructor.', ...
                  field);
        end
    end
end

% Create report object
o = class(o, 'report');
end