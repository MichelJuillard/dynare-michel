function b = uplus(a)

%@info:
%! @deftypefn {Function File} {@var{b} =} uplus (@var{a})
%! @anchor{@dynDate/uplus}
%! @sp 1
%! Overloads the uplus (unary addition) operator for the Dynare dates class (@ref{dynDate}). Increment the date by one year, quarter,
%! month or week depending on the frequency.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Dynare date object instantiated by @ref{dynDate}.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item b
%! Dynare date object instantiated by @ref{dynDate}.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @ref{dynDate}
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

if ~isa(a,'dynDate')
    error(['dynDate::uplus: Input argument ' inputname(1) ' must be a dynDate object.'])
end

b = dynDate(a);

switch b.freq
  case 1
    b.time(1) = b.time(1)+1;
  case 4
    if b.time(2)==4% Happy new year!
        b.time(1) = b.time(1)+1;
        b.time(2) = 1;
    else
        b.time(2) = b.time(2)+1;
    end
  case 12
    if b.time(2)==12% Happy new year!
        b.time(1) = b.time(1)+1;
        b.time(2) = 1;
    else
        b.time(2) = b.time(2)+1;
    end
  case 52
    if b.time(2)==52% Happy new year!
        b.time(1) = b.time(1)+1;
        b.time(2) = 1;
    else
        b.time(2) = b.time(2)+1;
    end
  otherwise
    error('dynDate::uplus: Unknown frequency!')
end

%@test:1
%$ addpath ../matlab
%$
%$ % Define some dates
%$ date_1 = '1950Q3';
%$ date_2 = '1950Q4';
%$ date_3 = '1950M3';
%$ date_4 = '1950M12';
%$ date_5 = '1950W3';
%$ date_6 = '1950W52';
%$ date_7 = 2000;
%$
%$ % Call the tested routine.
%$ d1 = dynDate(date_1);
%$ d2 = dynDate(date_2);
%$ d3 = dynDate(date_3);
%$ d4 = dynDate(date_4);
%$ d5 = dynDate(date_5);
%$ d6 = dynDate(date_6);
%$ d7 = dynDate(date_7);
%$ e1 = +d1;
%$ e2 = +d2;
%$ e3 = +d3;
%$ e4 = +d4;
%$ e5 = +d5;
%$ e6 = +d6;
%$ e7 = +d7;
%$
%$ % Check the results.
%$ t(1) = dyn_assert(e1.time,[1950 4]);
%$ t(2) = dyn_assert(e2.time,[1951 1]);
%$ t(3) = dyn_assert(e3.time,[1950 4]);
%$ t(4) = dyn_assert(e4.time,[1951 1]);
%$ t(5) = dyn_assert(e5.time,[1950 4]);
%$ t(6) = dyn_assert(e6.time,[1951 1]);
%$ t(7) = dyn_assert(e7.time,[2001 1]);
%$ T = all(t);
%@eof:1