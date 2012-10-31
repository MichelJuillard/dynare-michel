function c = max(a,b)

%@info:
%! @deftypefn {Function File} {@var{c} =} max (@var{a},@var{b})
%! @anchor{@dynDate/gt}
%! @sp 1
%! Overloads the max function for the Dynare dates class @ref{dynDate}.
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
%! Dynare date object instantiated by @ref{dynDate}. @var{c} is a copy of @var{a} if @var{a}>=@var{b}, a copy of @var{b} otherwise.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
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

verbose = 0;

if nargin~=2
    error('dynDate::min: I need exactly two input arguments!')
end

if ~( isa(a,'dynDate') && isa(b,'dynDate'))
    error(['dynDate::min: Input arguments ' inputname(1) 'and ' inputname(2) ' have to be a dynDate objects!'])
end

if verbose && a.freq~=b.freq
    error(['dynDate::min: Input arguments ' inputname(1) 'and ' inputname(2) ' have no common frequencies!'])
end

if a>=b
    c = a;
else
    c = b;
end

%@test:1
%$ addpath ../matlab
%$
%$ % Define some dates
%$ date_1 = 1950;
%$ date_2 = 2000;
%$ date_3 = '1950Q2';
%$ date_4 = '1950Q3';
%$ date_5 = '1950M1';
%$ date_6 = '1948M6';
%$
%$ % Call the tested routine.
%$ d1 = dynDate(date_1);
%$ d2 = dynDate(date_2);
%$ d3 = dynDate(date_3);
%$ d4 = dynDate(date_4);
%$ d5 = dynDate(date_5);
%$ d6 = dynDate(date_6);
%$ m1 = max(d1,d2);
%$ i1 = (m1==d2);
%$ m2 = max(d3,d4);
%$ i2 = (m2==d4);
%$ m3 = max(d5,d6);
%$ i3 = (m3==d5);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(i1,1);
%$ t(2) = dyn_assert(i2,1);
%$ t(3) = dyn_assert(i3,1);
%$ T = all(t);
%@eof:1
