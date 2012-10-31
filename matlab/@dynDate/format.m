function p = format(date)

%@info:
%! @deftypefn {Function File} {@var{p} =} format (@var{date})
%! @anchor{@dynDate/format}
%! @sp 1
%! Produces a formatted date from a Dynare date object instantiated by @ref{dynDate}.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item date
%! Dynare date object, instantiated by @ref{dynDate}.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item p
%! A string containing the formatted date (for instance, '2000', '2000Q3', '2000M9' or '2000W43').
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


if nargin~=1
    error('dynDate::format: I need exactly one input argument!')
end

if ~isa(date,'dynDate')
    error(['dynDate::format: Input argument ' inputname(1) ' has to be a dynDate object!'])
end

switch date.freq
  case 1
    p = num2str(date.time(1));
  case 4
    p = [num2str(date.time(1)) 'Q' num2str(date.time(2))];
  case 12
    p = [num2str(date.time(1)) 'M' num2str(date.time(2))];
  case 52
    p = [num2str(date.time(1)) 'W' num2str(date.time(2))];
  otherwise
    error('dynDate::format: Unkonwn frequency!')
end

%@test:1
%$ addpath ../matlab
%$
%$ % Define some dates
%$ date_1 = 1950;
%$ date_2 = '1950Q2';
%$ date_3 = '1950M10';
%$ date_4 = '1950W50';
%$
%$ % Call the tested routine.
%$ d1 = dynDate(date_1); DATE_1 = format(d1);
%$ d2 = dynDate(date_2); DATE_2 = format(d2);
%$ d3 = dynDate(date_3); DATE_3 = format(d3);
%$ d4 = dynDate(date_4); DATE_4 = format(d4);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(num2str(date_1),DATE_1);
%$ t(2) = dyn_assert(date_2,DATE_2);
%$ t(3) = dyn_assert(date_3,DATE_3);
%$ t(4) = dyn_assert(date_4,DATE_4);
%$ T = all(t);
%@eof:1
