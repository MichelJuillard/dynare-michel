function c = ge(a,b)

%@info:
%! @deftypefn {Function File} {@var{c} =} ge (@var{a},@var{b})
%! @anchor{@dynDate/ge}
%! @sp 1
%! Overloads the ge (greater or equal) operator for the Dynare dates class (@ref{dynDate}).
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
%! scalar integer equal to one if a>=b, 0 otherwise.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @ref{@@dynDate/gt}, @ref{@@dynDate/eq}
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

if a>b
    c=1;
else
    if a==b
        c=1;
    else
        c=0;
    end
end

%@test:1
%$ addpath ../matlab
%$
%$ % Define some dates
%$ date_1 = '1950Q3';
%$ date_2 = '1950Q3';
%$ date_3 = '1950Q1';
%$ date_4 = '1949Q2';
%$
%$ % Call the tested routine.
%$ d1 = dynDate(date_1);
%$ d2 = dynDate(date_2);
%$ d3 = dynDate(date_3);
%$ d4 = dynDate(date_4);
%$ i1 = (d1>=d2);
%$ i2 = (d3>=d4);
%$ i3 = (d4>=d2);
%$ i4 = (d1>=d4);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(i1,1);
%$ t(2) = dyn_assert(i2,1);
%$ t(3) = dyn_assert(i3,0);
%$ t(4) = dyn_assert(i4,1);
%$ T = all(t);
%@eof:1
