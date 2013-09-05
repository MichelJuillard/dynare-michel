function dd = pop(dd,a) % --*-- Unitary tests --*--

% pop method for dynDates class (removes a date)).

%@info:
%! @deftypefn {Function File} {@var{dd} =} pop (@var{dd}, @var{a})
%! @anchor{dynDates/pop}
%! @sp 1
%! Pop method for the Dynare dates class. Removes a date from a dynDates object, by default the laste date is removed.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item dd
%! Object instantiated by @ref{dynDates}.
%! @item a
%! Object instantiated by @ref{dynDate}, date to be removed.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dd
%! Object instantiated by @ref{dynDates}, without date (@var{a}).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2012 Dynare Team
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
    error(['dynDates::pop: Input argument ' inputname(dd) ' has to be a dynDates object.'])
end

if nargin<2
    % Remove last date
    dd.time = dd.time(1:end-1,:);
    dd.ndat = dd.ndat-1;
    return
end

if nargin>1 && ~(isa(a,'dynDate') || ischar(a))
    error(['dynDates::pop: Input argument ' inputname(a) ' has to be a dynDate object or a string (formatted date).'])
end

if ~isa(a,'dynDate')
    a = dynDate(a);
end

idx = find(~all(transpose(bsxfun(@eq,dd.time,a.time))));    
dd.time = dd.time(idx,:);
dd.ndat = size(dd.time,1);

%@test:1
%$ % Define some dates
%$ B1 = '1953Q4';
%$ B2 = '1950Q2';
%$ B3 = '1950Q1';
%$ B4 = '1945Q3';
%$ B5 = '2009Q2';
%$
%$ % Define expected results
%$ e.time = [1945 3; 1950 1; 1950 2; 1953 4; 2009 2];
%$ e.freq = 4;
%$ e.ndat = 5;
%$
%$ % Call the tested routine
%$ d = dynDates(B4,B3,B2,B1);
%$ d = d.append(dynDate(B5));
%$ f = d.pop();
%$ t(1) = dyn_assert(f.time,e.time(1:end-1,:));
%$ t(2) = dyn_assert(f.freq,e.freq);
%$ t(3) = dyn_assert(f.ndat,e.ndat-1);
%$ f = d.pop(B1);
%$ t(4) = dyn_assert(f.time,[1945 3; 1950 1; 1950 2; 2009 2]);
%$ t(5) = dyn_assert(f.freq,e.freq);
%$ t(6) = dyn_assert(f.ndat,e.ndat-1);
%$ f = d.pop(dynDate(B1));
%$ t(7) = dyn_assert(f.time,[1945 3; 1950 1; 1950 2; 2009 2]);
%$ t(8) = dyn_assert(f.freq,e.freq);
%$ t(9) = dyn_assert(f.ndat,e.ndat-1);
%$
%$ % Check the results.
%$ T = all(t);
%@eof:1