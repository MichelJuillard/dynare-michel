function p = page(varargin)
%function p = page(varargin)
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

p = struct;
p.sections = sections();

switch nargin
    case 0
        p = class(p, 'page');
    case 1
        assert(isa(varargin{1}, 'page'), ['Page constructor: the only valid ' ...
            'arguments are page objects']);
        p = varargin{1};
    otherwise
        error('Page constructor: invalid number of arguments');
end
end