function sp = colon(a,b)

%@info:
%! @deftypefn {Function File} {@var{sp} =} colon (@var{a},@var{b})
%! @anchor{@dynDate/colon}
%! @sp 1
%! Overloads the colon operator for the Dynare Dates class (@ref{dynDate}). Creates a @ref{dynTime} object.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Dynare date object instantiated by @ref{dynDate}, initial date.
%! @item b
%! Dynare date object instantiated by @ref{dynDate}, last date.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item c
%! Dynare Time object instantiated by @ref{dynTime}.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @ref{dynTime}, @ref{@@dynTime/setFreq}, @ref{@@dynTime/setTime}, @ref{@@dynTime/setSize}
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

if nargin~=2
    error('dynTime::colon: I need exactly two input arguments!')
end

if ~( isa(a,'dynDate') && isa(b,'dynDate'))
    error(['dynTime::colon: Input arguments ' inputname(1) 'and ' inputname(2) ' have to be a dynDate objects!'])
end

if a.freq~=b.freq
    error(['dynTime::colon: Input arguments ' inputname(1) 'and ' inputname(2) ' must have common frequency!'])
end

if a>b
    error(['dynTime::colon: ' inputname(1) ' must precede ' inputname(2) '!' ])
end

if a==b% Time range with only one date.
    sp = dynTime();
    sp = sp.setFreq(a.freq);
    sp = sp.setSize(n+1);
    sp = sp.setTime(1,a.time);
else
    n = b-a;
    sp = dynTime();
    sp = sp.setFreq(a.freq);
    sp = sp.setSize(n+1);
    sp = sp.setTime(1,a.time);
    for t=2:n+1
        a = +a;
        sp = sp.setTime(t,a.time);
    end
end

%@test:1
%$ addpath ../matlab
%$
%$ % Define two dates
%$ date_1 = '1950Q2';
%$ date_2 = '1951Q4';
%$
%$ % Define expected results.
%$ e.freq = 4;
%$ e.time = [1950 2; 1950 3; 1950 4; 1951 1; 1951 2; 1951 3; 1951 4];
%$
%$ % Call the tested routine.
%$ d1 = dynDate(date_1);
%$ d2 = dynDate(date_2);
%$ d3 = d1:d2;
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d3.time,e.time);
%$ t(2) = dyn_assert(d3.freq,e.freq);
%$ T = all(t);
%@eof:1
