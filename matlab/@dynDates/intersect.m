function C = intersect(A,B)

%@info:
%! @deftypefn {Function File} {@var{C} =} intersect (@var{A},@var{B})
%! @anchor{@dynDates/intersect}
%! @sp 1
%! C of B and A.
%! if A and B are not disjoints.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! @ref{dynDates} object.
%! @item B
%! @ref{dynDates} object.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item C
%! @ref{dynDates} object.
%! @end table
%! @end deftypefn
%@eod:

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

if ~isa(A,'dynDates') || ~isa(B,'dynDates')
    error(['dynDates::plus: Input arguments ''' inputname(1) ''' and ''' inputname(2) ''' must be dynDates objects!'])
end

if eq(A,B)
    C = A;
    return
end

if ~isequal(A.freq,B.freq)
    C = dynDates();
    return
end

time = intersect(A.time,B.time,'rows');

C = dynDates();
if isempty(time)
    return
end

C.freq = A.freq;
C.time = time;
C.ndat = rows(time); 

%@test:1
%$ % Define some dynDates objects
%$ d1 = dynDate('1950Q1'):dynDate('1969Q4') ;
%$ d2 = dynDate('1960Q1'):dynDate('1969Q4') ;
%$ d3 = dynDate('1970Q1'):dynDate('1979Q4') ;
%$
%$ % Call the tested routine.
%$ c1 = intersect(d1,d2);
%$ c2 = intersect(d1,d3);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(c1==d2,1);
%$ t(2) = dyn_assert(isempty(c2),1);
%$ T = all(t);
%@eof:1