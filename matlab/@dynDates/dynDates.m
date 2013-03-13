function dd = dynDates(varargin)

%@info:
%! @deftypefn {Function File} {@var{dd} =} dynDates (@var{a},@var{b},...)
%! @anchor{dynDates}
%! @sp 1
%! Constructor for the Dynare dates class (unordered sequence of dates).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! String, date.
%! @item b
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dd
%! Dynare dates object.
%! @end table
%! @sp 1
%! @strong{Properties}
%! @sp 1
%! The constructor defines the following properties:
%! @sp 1
%! @table @ @var
%! @item ndate
%! Scalar integer, the number of dates.
%! @item freq
%! Scalar integer, the frequency of the time series. @var{freq} is equal to 1 if data are on a yearly basis or if
%! frequency is unspecified. @var{freq} is equal to 4 if data are on a quaterly basis. @var{freq} is equal to
%! 12 if data are on a monthly basis. @var{freq} is equal to 52 if data are on a weekly basis.
%! @item time
%! Array of integers (nobs*2). The first column defines the years associated to each date. The second column,
%! depending on the frequency, indicates the week, month or quarter numbers. For yearly data or unspecified frequency
%! the second column is filled by ones.
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

dd = struct;

dd.ndat = 0;
dd.freq = [];
dd.time = [];

dd = class(dd,'dynDates');

switch nargin
  case 0
    % Returns an empty object
    return
  case 1
    if isa(varargin{1},'dynDates')
        % Returns a copy of the input argument
        dd = varargin{1};
    elseif ischar(varargin{1}) || isa(varargin{1},'dynDate')
        tmp = dynDate(varargin{1});
        dd.ndat = 1;
        dd.freq = tmp.freq;
        dd.time = tmp.time;
    else
        error('dynDates:: Wrong calling sequence of the constructor!')
    end
  otherwise
    tmp = dynDate(varargin{1});
    dd.ndat = nargin;
    dd.time = NaN(dd.ndat,2);
    dd.freq = tmp.freq;
    dd.time(1,:) = tmp.time;
    for i=2:dd.ndat
        tmp = dynDate(varargin{i});
        if ~isequal(dd.freq,tmp.freq)
            error(['dynDates:: The frequency declared in input argument number ' int2str(i) ' is different from the frequency declared in the first input argument!'])
        end
        dd.time(i,:) = tmp.time;
    end
end

%@test:1
%$ % Define some dates
%$ B1 = '1945Q3';
%$ B2 = '1950Q2';
%$ B3 = '1950Q1';
%$ B4 = '1953Q4';
%$
%$ % Define expected results.
%$ e.time = [1945 3; 1950 2; 1950 1; 1953 4];
%$ e.freq = 4;
%$ e.ndat = 4;
%$
%$ % Call the tested routine.
%$ d = dynDates(B1,B2,B3,B4);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d.time,e.time);
%$ t(2) = dyn_assert(d.freq,e.freq);
%$ t(3) = dyn_assert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define some dates
%$ B1 = '1945M3';
%$ B2 = '1950M2';
%$ B3 = '1950M10';
%$ B4 = '1953M12';
%$
%$ % Define expected results.
%$ e.time = [1945 3; 1950 2; 1950 10; 1953 12];
%$ e.freq = 12;
%$ e.ndat = 4;
%$
%$ % Call the tested routine.
%$ d = dynDates(B1,B2,B3,B4);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d.time,e.time);
%$ t(2) = dyn_assert(d.freq,e.freq);
%$ t(3) = dyn_assert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define some dates
%$ B1 = '1945';
%$ B2 = '1950';
%$ B3 = '1950';
%$ B4 = '1953';
%$
%$ % Define expected results.
%$ e.time = [1945 1; 1950 1; 1950 1; 1953 1];
%$ e.freq = 1;
%$ e.ndat = 4;
%$
%$ % Call the tested routine.
%$ d = dynDates(B1,B2,B3,B4);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d.time,e.time);
%$ t(2) = dyn_assert(d.freq,e.freq);
%$ t(3) = dyn_assert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define some dates
%$ B1 = '1945Q1';
%$ B2 = '1950Q3';
%$ B3 = '1950M10';
%$ B4 = '1953Q1';
%$
%$
%$ % Call the tested routine.
%$ try
%$     d = dynDates(B1,B2,B3,B4);
%$     t(1) = 0;
%$     T = 0;
%$ catch
%$     % Expected issue...
%$     t(1) = 1;
%$     T = 1;
%$ end
%@eof:4