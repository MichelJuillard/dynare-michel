function dd = append(dd,a)
% append method for dynDates class.

%@info:
%! @deftypefn {Function File} {@var{dd} =} sort (@var{dd}, @var{a})
%! @anchor{dynDates/append}
%! @sp 1
%! Append method for the Dynare dates class.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item dd
%! Object instantiated by @ref{dynDates}.
%! @item a
%! Object instantiated by @ref{dynDate}.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dd
%! Object instantiated by @ref{dynDates}, with an additional date (@var{a}).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2012-2013 Dynare Team
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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

if ~isa(dd,'dynDates')
    error(['dynDates::append: Input argument ' inputname(dd) ' has to be a dynDates object.'])
end

if ~(isa(a,'dynDate') || isa(a,'dynDates') || ischar(a))
    error(['dynDates::append: Input argument ' inputname(a) ' has to be ' ...
                        'a dynDate object or a dynDates object or a string (formatted date).'])
end

if isempty(a)
    return
end

if isa(a,'dynDate')
    dd.time = [dd.time; a.time];
    dd.ndat = dd.ndat+1;
elseif isa(a,'dynDates')
    dd.time = [dd.time; a.time];
    dd.ndat = dd.ndat+a.ndat;
else
    tmp = dynDate();
    dd.time = [dd.time; tmp(a).time];
    dd.ndat = dd.ndat+1;
end

%@test:1
%$ % Define some dates
%$ B1 = '1953Q4';
%$ B2 = '1950Q2';
%$ B3 = '1950Q1';
%$ B4 = '1945Q3';
%$ B5 = '2009Q2';
%$
%$ % Define expected results.
%$ e.time = [1945 3; 1950 1; 1950 2; 1953 4; 2009 2];
%$ e.freq = 4;
%$ e.ndat = 5;
%$
%$ % Call the tested routine.
%$ d = dynDates(B4,B3,B2,B1);
%$ d = d.append(dynDate(B5));
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d.time,e.time);
%$ t(2) = dyn_assert(d.freq,e.freq);
%$ t(3) = dyn_assert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define some dates
%$ B1 = '1953Q4';
%$ B2 = '1950Q2';
%$ B3 = '1950Q1';
%$ B4 = '1945Q3';
%$ B5 = '2009q2';
%$
%$ % Define expected results.
%$ e.time = [1945 3; 1950 1; 1950 2; 1953 4; 2009 2];
%$ e.freq = 4;
%$ e.ndat = 5;
%$
%$ % Call the tested routine.
%$ d = dynDates(B4,B3,B2,B1);
%$ d = d.append(B5);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d.time,e.time);
%$ t(2) = dyn_assert(d.freq,e.freq);
%$ t(3) = dyn_assert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:2