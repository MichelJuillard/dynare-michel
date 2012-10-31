function c = eq(a,b)

%@info:
%! @deftypefn {Function File} {@var{c} =} eq (@var{a},@var{b})
%! @anchor{@dynDate/eq}
%! @sp 1
%! Overloads the eq (equal) operator for the Dynare dates class (@ref{dynDate}).
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
%! scalar integer equal to one if a==b, 0 otherwise.
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
    disp(['dynDate::eq: Input arguments ' inputname(1) 'and ' inputname(2) ' have no common frequencies!'])
    c = 0;
end

c = isequal(a.time,b.time);

%@test:1
%$ addpath ../matlab
%$
%$ % Define some dates
%$ date_1 = 1950;
%$ date_2 = '1950Q2';
%$ date_3 = '1950M10';
%$ date_4 = '1950W50';
%$ date_5 = '1950W32';
%$
%$ % Call the tested routine.
%$ d1 = dynDate(date_1);
%$ d2 = dynDate(date_2);
%$ d3 = dynDate(date_3);
%$ d4 = dynDate(date_4);
%$ d5 = dynDate(date_5);
%$ i1 = (d1==d2);
%$ i2 = (d2==d2);
%$ i3 = (d4==d5);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(i1,0);
%$ t(2) = dyn_assert(i2,1);
%$ t(3) = dyn_assert(i3,0);
%$ T = all(t);
%@eof:1

