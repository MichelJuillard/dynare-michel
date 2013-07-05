% --*-- Unitary tests --*--
function C = plus(A,B)

%@info:
%! @deftypefn {Function File} {@var{C} =} plus (@var{A},@var{B})
%! @anchor{@dynDates/plus}
%! @sp 1
%! Overloads the plus (addition) operator for the @ref{dynDates} class. Combines two dynDates objects, A and B, without removing repetitions
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

if isempty(B)
    C = A;
    return
end

if isempty(A)
    C = B;
    return
end

if ~isequal(A.freq,B.freq)
    error(['dynDates::plus: Input arguments ''' inputname(1) ''' and ''' inputname(2) ''' must have common frequencies!'])
end

C = dynDates();

C.freq = A.freq;
C.time = [A.time; B.time];
C.ndat = A.ndat+B.ndat;

%@test:1
%$ % Define some dynDates objects
%$ d1 = dynDate('1950Q1'):dynDate('1959Q4') ;
%$ d2 = dynDate('1960Q1'):dynDate('1969Q4') ;
%$ d3 = dynDate('1970Q1'):dynDate('1979Q4') ;
%$
%$ % Call the tested routine.
%$ e1 = d1+d2;
%$ e2 = d1+d2+d3;
%$
%$ % Expected results.
%$ f1 = dynDate('1950Q1'):dynDate('1969Q4');
%$ f2 = dynDate('1950Q1'):dynDate('1979Q4');
%$
%$ % Check the results.
%$ t(1) = dyn_assert(e1==f1,1);
%$ t(2) = dyn_assert(e2==f2,1);
%$ T = all(t);
%@eof:1