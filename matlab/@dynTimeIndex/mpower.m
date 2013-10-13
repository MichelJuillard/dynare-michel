function C = mpower(A,B) % --*-- Unitary tests --*--

% C = mpower(A,B)
%
% Overloads binary mpower operator (^).
%
% INPUTS :
%  * A, dynTimeIndex object.
%  * B, integer scalar.
%
% OUTPUTS :
%  * C, dynTimeIndex object.
%
% EXAMPLE 1 :
%
%  >> B = dynTimeIndex()-1;
%  >> B
%  B = <dynTimeIndex: -1>
%  >> B^4
%  ans = <dynTimeIndex: -4>
%  >>
%
%  EXAMPLE 2 :
%  This method can be used to apply the lead and lag methods an arbitrary number of times to a dynSeries object. For instance, if 
%  ts is a dynSeries object, and if we define
% 
%  >> B = dynTimeIndex()-1;
%  >> F = dynTimeIndex()+1;
%
%  B and F can be used as lag and lead operators and the following syntax:
%
%  >> us = ts(F^2);
%
%  is equivalent to
%
%  >> us = ts.lead(2)
%
%  or
%
%  >> us = ts.lead.lead
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
    error(['dynTimeIndex::mpower: Input arguments ''' inputname(1) ''' and ''' inputname(2) ''' must be a dynTimeIndex object and an integer!'])
end

C = struct();
C.index = A.index*B;
C = class(C,'dynTimeIndex');

%@test:1
%$ a = dynTimeIndex()+1;
%$ b = a^2;
%$ t(1) = isa(b,'dynTimeIndex');
%$ t(2) = isequal(b.index,int8(2));
%$ T = all(t);
%@eof:1