function oa = addObj(oa, varargin)
%function oa = addObj(oa, varargin)

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

assert(nargin >= 2 && nargin <= 3)
assert(isa(oa, 'objArray'), 'First argument must be an objArray');
assert(isobject(varargin{1}), 'Optional 2nd arg must be an object');
if nargin == 3
    assert(isnumeric(varargin{2}), 'Optional 3rd arg must be an index');
end

if nargin == 2
    oa.objs{end+1} = varargin{1};
elseif nargin == 3
    oa.objs{varargin{2}} = varargin{1};
end