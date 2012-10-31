function c = lt(a,b)

%@info:
%! @deftypefn {Function File} {@var{c} =} lt (@var{a},@var{b})
%! @anchor{@dynDate/lt}
%! @sp 1
%! Overloads the lt (less than) operator for the Dynare dates class (@ref{dynDate}).
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
%! scalar integer equal to one if a<b, 0 otherwise.
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
    error('dynDate::eq: I need exactly two input arguments!')
end

if ~( isa(a,'dynDate') && isa(b,'dynDate'))
    error(['dynDate::eq: Input arguments ' inputname(1) 'and ' inputname(2) ' have to be a dynDate objects!'])
end

if verbose && a.freq~=b.freq
    error(['dynDate::eq: Input arguments ' inputname(1) 'and ' inputname(2) ' have no common frequencies!'])
end

if a.time(1)<b.time(1)
    c = 1;
elseif isequal(a.time(1),b.time(1))
    if a.time(2)<b.time(2)
        c = 1;
    else
        c = 0;
    end
else
    c = 0;
end

%@test:1
%$ addpath ../matlab
%$
%$ % Define some dates
%$ date_1 = 1950;
%$ date_2 = '1950Q2';
%$ date_3 = '1950Q3';
%$ date_4 = '1950Q1';
%$ date_5 = '1949Q2';
%$
%$ % Call the tested routine.
%$ d1 = dynDate(date_1);
%$ d2 = dynDate(date_2);
%$ d3 = dynDate(date_3);
%$ d4 = dynDate(date_4);
%$ d5 = dynDate(date_5);
%$ i1 = (d2<d3);
%$ i2 = (d3<d4);
%$ i3 = (d4<d2);
%$ i4 = (d5<d4);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(i1,1);
%$ t(2) = dyn_assert(i2,0);
%$ t(3) = dyn_assert(i3,1);
%$ t(4) = dyn_assert(i4,1);
%$ T = all(t);
%@eof:1
