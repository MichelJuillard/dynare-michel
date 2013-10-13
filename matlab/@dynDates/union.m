function C = union(A,varargin) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{C} =} union (@var{A},@var{B})
%! @anchor{dynDates/union}
%! @sp 1
%! Union method for the Dynare dates class (removes repetitions if any). Dates in C are sorted in increasing order.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Object instantiated by @ref{dynDates}.
%! @item B
%! Object instantiated by @ref{dynDates}.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item C
%! Object instantiated by @ref{dynDates}.
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

if ~isa(A,'dynDates')  
    error(['dynDates::union: Input argument ''' inputname(1) '''   has to be a dynDates object!'])
end

C = A;

if ~length(varargin)
    return;
end

for i=1:length(varargin)
    if isa(varargin{i},'dynDates')
        C = C + varargin{i};
    else
        error(['dynDates::union: Input argument ''' inputname(i) '''   has to be a dynDates object!'])
    end
end

C = sort(unique(C));

%@test:1
%$ % Define some dynDates objects
%$ d1 = dynDate('1950Q1'):dynDate('1959Q4') ;
%$ d2 = dynDate('1960Q1'):dynDate('1969Q4') ;
%$ d3 = dynDate('1970Q1'):dynDate('1979Q4') ;
%$
%$ % Call the tested routine.
%$ e1 = union(d1);
%$ e2 = union(d1,d2);
%$ e3 = union(d1,d2,d3);
%$ e4 = union(d1,d2,d3,d2+d3);
%$ e5 = union(d1,d2,d3,d2);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(e1==d1,1);
%$ t(2) = dyn_assert(e2==d1+d2,1);
%$ t(3) = dyn_assert(e3==d1+d2+d3,1);
%$ t(4) = dyn_assert(e4==d1+d2+d3,1);
%$ t(5) = dyn_assert(e5==d1+d2+d3,1);
%$ T = all(t);
%@eof:1