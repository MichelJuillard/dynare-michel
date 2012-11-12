function ts = dynSeries(a,b,c,d)

%@info:
%! @deftypefn {Function File} {@var{ts} =} dynSeries (@var{a},@var{b},@var{c},@var{d})
%! @anchor{dynSeries}
%! @sp 1
%! Constructor for the Dynare time series class.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! T*1 vector or T*N matrix of data.
%! @item b
%! Initial date. For Quaterly, Monthly or Weekly data, b must be a string. For yearly data or if the frequence is not
%! defined b must be an integer.
%! @item c
%! N*q array of characters. Names of the N time series.
%! @item d
%! N*p array of characters. TeX names of the N time series.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item ts
%! Dynare time series object.
%! @end table
%! @sp 1
%! @strong{Properties}
%! @sp 1
%! The constructor defines the following properties:
%! @sp 1
%! @table @ @var
%! @item data
%! Array of doubles (nobs*vobs).
%! @item nobs
%! Scalar integer, the number of observations.
%! @item vobs
%! Scalar integer, the number of variables.
%! @item name
%! Array of chars (nvobs*n), names of the variables.
%! @item tex
%! Array of chars (nvobs*n), tex names of the variables.
%! @item freq
%! Scalar integer, the frequency of the time series. @var{freq} is equal to 1 if data are on a yearly basis or if
%! frequency is unspecified. @var{freq} is equal to 4 if data are on a quaterly basis. @var{freq} is equal to
%! 12 if data are on a monthly basis. @var{freq} is equal to 52 if data are on a weekly basis.
%! @item time
%! Array of integers (nobs*2). The first column defines the years associated to each observation. The second column,
%! depending on the frequency, indicates the week, month or quarter numbers. For yearly data or unspecified frequency
%! the second column is filled by ones.
%! @item init
%! Row vector of integers (1*2) indicating the year and the week, month or quarter of the first observation. @var{init}
%! is the first row of @var{time}.
%! @item last
%! Row vector of integers (1*2) indicating the year and the week, month or quarter of the last observation. @var{init}
%! is the first row of @var{time}.
%! @end table
%! @sp 1
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @ref{@@dynTime/dynTime}, @ref{@@dynTime/setTime}, @ref{@@dynTime/setFreq} 
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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

ts = struct;

ts.data = [];
ts.nobs = 0;
ts.vobs = 0;
ts.name = [];
ts.tex  = [];
ts.freq = [];
ts.Time = dynTime();

ts = class(ts,'dynSeries');

switch nargin
  case 0
    return
  case {2,4}
    if nargin==2
        c = [];
        d = [];
    end
    % Get data, number of observations and number of variables.
    ts.data = a;
    ts.nobs = size(a,1);
    ts.vobs = size(a,2);
    % Get the first date and set the frequency.
    if ~isempty(b)
        if ischar(b)% Weekly, Monthly or Quaterly data.
            quaterly = findstr('Q',b);
            monthly  = findstr('M',b);
            weekly   = findstr('W',b);
            if ~isempty(quaterly)
                ts.freq = 4;
            end
            if ~isempty(monthly)
                ts.freq = 12;
            end
            if ~isempty(weekly)
                ts.freq = 52;
            end
            if isempty(quaterly) && isempty(monthly) && isempty(weekly)
                error('dynSeries:: Using a string as a second input argument, I can only handle weekly (W), monthly (M) or quaterly (Q) data!');
            end
        else% If b is not a string then yearly data are assumed.
            ts.freq = 1;
        end
        ts.Time = ts.Time.setFreq(ts.freq);
        ts.Time = ts.Time.setTime(dynDate(b):dynDate(b)+ts.nobs);
    else% If b is empty.
        ts.freq = 1;
        ts.Time = ts.Time.setFreq(1);
        ts.Time = ts.Time.setTime([transpose(1:ts.nobs) ones(ts.nobs,1)]);
    end
    % Get the names of the variables.
    if ~isempty(c)
        if ts.vobs==size(c,1)
            ts.name = c;
        else
            error('dynSeries:: The number of declared names does not match the number of variables!')
        end
    else
        for i=1:ts.vobs
            ts.name = char(ts.name,'--NA--');
        end
    end
    if ~isempty(d)
        if ts.vobs==size(d,1)
            ts.tex = d;
        else
            error('dynSeries:: The number of declared tex names does not match the number of variables!')
        end
    else
        for i=1:ts.vobs
            ts.tex = char(ts.tex,'--NA--');
        end
    end
  otherwise
    error('dynSeries:: Can''t instantiate the class, wrong calling sequence!')
end


%@test:1
%$ % Test if we can instantiate an empty dynSeries object.
%$ try
%$     ts = dynSeries();
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$ disp('coucou')
%$ T = all(t);
%@eof:1

%@test:2
%$ addpath ../matlab
%$ % Define a data set.
%$ A = transpose(1:10);
%$
%$ % Define initial date
%$ B1 = 1950;
%$ B2 = '1950Q2';
%$ B3 = '1950M10';
%$ B4 = '1950W50';
%$
%$ % Define expected results.
%$ e1.Time = transpose([1950 1951 1952 1953 1954 1955 1956 1957 1958 1959]);
%$ e1.freq = 1;
%$ e2.Time = char('1950Q2','1950Q3','1950Q4','1951Q1','1951Q2','1951Q3','1951Q4','1952Q1','1952Q2','1952Q3');
%$ e2.freq = 4;
%$ e3.Time = char('1950M10','1950M11','1950M12','1951M1','1951M2','1951M3','1951M4','1951M5','1951M6','1951M7');
%$ e3.freq = 12;
%$ e4.Time = char('1950W50','1950W51','1950W52','1951W1','1951W2','1951W3','1951W4','1951W5','1951W6','1951W7');
%$ e4.freq = 52;
%$
%$ % Call the tested routine.
%$ ts1 = dynSeries(A,B1);
%$ ts2 = dynSeries(A,B2);
%$ ts3 = dynSeries(A,B3);
%$ ts4 = dynSeries(A,B4);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(getTime(ts1),e1.Time);
%$ t(2) = dyn_assert(getTime(ts2),e2.Time);
%$ t(3) = dyn_assert(getTime(ts3),e3.Time);
%$ t(4) = dyn_assert(getTime(ts4),e4.Time);
%$ t(5) = dyn_assert(ts1.freq,e1.freq);
%$ t(6) = dyn_assert(ts2.freq,e2.freq);
%$ t(7) = dyn_assert(ts3.freq,e3.freq);
%$ t(8) = dyn_assert(ts4.freq,e4.freq);
%$ T = all(t);
%@eof:1