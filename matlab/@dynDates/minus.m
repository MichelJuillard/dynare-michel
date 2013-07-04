% --*-- Unitary tests --*--
function C = minus(A,B)

%@info:
%! @deftypefn {Function File} {@var{C} =} minus (@var{A},@var{B})
%! @anchor{@dynDates/minus}
%! @sp 1
%! Overloads the minus (soustraction) operator for the @ref{dynDates} class. C is the relative complement of B in A.
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

if isempty(B)
    C = A;
    return
end

if isempty(A)
    C = dynDates();
    return
end

if ~isequal(A.freq,B.freq)
    C = A;
    return
end

D = intersect(A,B);

if isempty(D)
    C = A;
else
    C = dynDates();
    C.freq = A.freq;
    C.time = setdiff(A.time,D.time,'rows');
    C.ndat = rows(C.time);
end

%@test:1
%$ % Define some dynDates objects
%$ d1 = dynDate('1950Q1'):dynDate('1959Q4') ;
%$ d2 = dynDate('1960Q1'):dynDate('1969Q4') ;
%$ d3 = d1+d2;
%$
%$ % Call the tested routine.
%$ e1 = d1-d2;
%$ e2 = d3-d1;
%$
%$ % Check the results.
%$ t(1) = dyn_assert(e1==d1,1);
%$ t(2) = dyn_assert(e2==d2,1);
%$ T = all(t);
%@eof:1