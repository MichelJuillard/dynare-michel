function C = plus(A,B) % --*-- Unitary tests --*--

% C = plus(A,B)
%
% Overloads binary plus operator.
%
% INPUTS:
%  * A, dynTimeIndex object.
%  * B, integer scalar.
%
% OUTPUTS:
%  * C, dynTimeIndex object.
%
% EXAMPLE:
%
%  >> t = dynTimeIndex();
%  >> t.index
%
%  ans =
%
%      0
%
%  >> s = t+1;
%  >> s.index
%
%  ans =
%
%      1
%

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

if ~(isa(A,'dynTimeIndex') || isint(B))
    error(['dynTimeIndex::plus: Input arguments ''' inputname(1) ''' and ''' inputname(2) ''' must be a dynTimeIndex object and an integer!'])
end

C = struct();
C.index = A.index+B;
C = class(C,'dynTimeIndex');

%@test:1
%$ a = dynTimeIndex();
%$ b = a+1;
%$ t(1) = isa(b,'dynTimeIndex');
%$ t(2) = isequal(b.index,int8(1));
%$ T = all(t);
%@eof:1