function c = plus(a,b)

%@info:
%! @deftypefn {Function File} {@var{c} =} plus (@var{a},@var{b})
%! @anchor{@dynDates/plus}
%! @sp 1
%! Overloads the plus (addition) operator for the Dynare dates class (@ref{dynDates}). Given an initial date @var{a},
%! computes a new date @var{c} by adding the number of periods @var{b}.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Dynare date object instantiated by @ref{dynDates}.
%! @item b
%! Positive scalar integer, the number of periods
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item c
%! Dynare date object instantiated by @ref{dynDates}.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @ref{@@dynDates/eq}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011 Dynare Team
% stephane DOT adjemian AT univ DASH lemans DOT fr
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

if ~isa(a,'dynDates')
    error(['dynDates::plus: Input argument ' inputname(1) ' must be a dynDates object!'])
end

if b<0
    error(['dynDates::plus: Input argument ' inputname(2) ' must be a positive integer'])
end


if b==0
    c = a;
    return
end

switch a.freq
  case 1
    c = a.time(1)-b.time(1);
  case {4,12,52}
    c = a.time(2)-b.time(2) + (a.time(1)-b.time(1))*a.freq;
  otherwise
    error('dynDates::minus: Unknown frequency!')
end

%@test:1
%$ addpath ../matlab
%$
%$ % Define some dates
%$ date_0_1 = 1950;
%$ date_0_2 = 1950;
%$ date_0_3 = 1940;
%$ date_1_1 = '1950Q4';
%$ date_1_2 = '1950Q1';
%$ date_1_3 = '1940Q3';
%$ date_2_1 = '2000M3';
%$ date_2_2 = '1998M8';
%$
%$ % Call the tested routine.
%$ d_0_1 = dynDates(date_0_1);
%$ d_0_2 = dynDates(date_0_2);
%$ d_0_3 = dynDates(date_0_3);
%$ d_1_1 = dynDates(date_1_1);
%$ d_1_2 = dynDates(date_1_2);
%$ d_1_3 = dynDates(date_1_3);
%$ d_2_1 = dynDates(date_2_1);
%$ d_2_2 = dynDates(date_2_2);
%$ e1 = d_0_1-d_0_2;
%$ e2 = d_0_1-d_0_3;
%$ e3 = d_1_1-d_1_2;
%$ e4 = d_1_1-d_1_3;
%$ e5 = d_2_1-d_2_2;
%$
%$ % Check the results.
%$ t(1) = dyn_assert(e1,0);
%$ t(2) = dyn_assert(e2,10);
%$ t(3) = dyn_assert(e3,3);
%$ t(4) = dyn_assert(e4,41);
%$ t(4) = dyn_assert(e5,19);
%$ T = all(t);
%@eof:1