function e = elements(varargin)
%function e = elements(varargin)
% Elements Class Constructor
%
% INPUTS
%   Optional elements object
%
% OUTPUTS
%   elements object
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

switch nargin
    case 0
        e = class(struct, 'elements', objArray());
    case 1
        assert(isa(varargin{1}, 'elements'), ...
            ['With one arg to elements constructor, you must pass an ' ...
            'elements object or a char.']);
        e = varargin{1};
    otherwise
        error('Elements constructor: invalid number of arguments');
end
end