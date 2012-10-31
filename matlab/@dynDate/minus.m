function c = minus(a,b)

%@info:
%! @deftypefn {Function File} {@var{c} =} minus (@var{a},@var{b})
%! @anchor{@dynDate/minus}
%! @sp 1
%! Overloads the minus (soustraction) operator for the Dynare dates class (@ref{dynDate}). Depending on the frequency, computes the number
%! of years, quarters, months, weeks between two dates @var{a} and @var{b} (it is assumed that @var{a}>@var{B}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Dynare date object instantiated by @ref{dynDate}.
%! @item b
%! Dynare date object instantiated by @ref{dynDate}.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item c
%! Integer scalar, the number of years, quarters, months or weeks between @var{a} and @var{B}.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @ref{@@dynDate/eq},@ref{@@dynDate/lt}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011 Dynare Team
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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

if ~( isa(a,'dynDate') && isa(b,'dynDate') )
    error(['dynDate::minus: Input arguments ' inputname(1) ' and ' inputname(2) ' must be dynDate objects!'])
end

if a.freq~=b.freq
    error(['dynDate::minus: ' inputname(1) ' and ' inputname(2) ' must have common frequency!'])
end

if a<b
    error(['dynDate::minus: ' inputname(1) ' must be posterior to ' inputname(2) '!'])
end

if a==b
    c = 0;
    return
end

switch a.freq
  case 1
    c = a.time(1)-b.time(1);
  case {4,12,52}
    c = a.time(2)-b.time(2) + (a.time(1)-b.time(1))*a.freq;
  otherwise
    error('dynDate::minus: Unknown frequency!')
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
%$ d_0_1 = dynDate(date_0_1);
%$ d_0_2 = dynDate(date_0_2);
%$ d_0_3 = dynDate(date_0_3);
%$ d_1_1 = dynDate(date_1_1);
%$ d_1_2 = dynDate(date_1_2);
%$ d_1_3 = dynDate(date_1_3);
%$ d_2_1 = dynDate(date_2_1);
%$ d_2_2 = dynDate(date_2_2);
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