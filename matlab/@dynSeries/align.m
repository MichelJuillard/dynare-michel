function [a,b] = align(a, b) % --*-- Unitary tests --*--
    
%@info:
%! @deftypefn {Function File} {[@var{a}, @var{b}] =} align (@var{a}, @var{b})
%! @anchor{dynSeries/align}
%! @sp 1
%! If dynSeries objects @var{a} and @var{b} are defined on different time ranges, extend @var{a} and/or
%! @var{b} with NaNs so that they are defined on the same time range.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Object instantiated by @ref{dynSeries}.
%! @item b
%! Object instantiated by @ref{dynSeries}.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Object instantiated by @ref{dynSeries}.
%! @item b
%! Object instantiated by @ref{dynSeries}.
%! @end table
%! @end deftypefn
%@eod:

% Copyright (C) 2013 Dynare Team
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

if ~isequal(a.freq,b.freq)
    error(['dynSeries::merge_dates: ''' inputname(1) ''' and ''' inputname(2) ''' dynSeries objects must have common frequencies!'])
end

init = min(a.init,b.init);

time_range_of_a = a.init:a.init+a.nobs;
time_range_of_b = b.init:b.init+b.nobs;

last_a = time_range_of_a(a.nobs);
last_b = time_range_of_b(b.nobs);

common_time_range = intersect(time_range_of_a,time_range_of_b);

if isempty(common_time_range)
    error(['dynSeries::merge_dates: ''' inputname(1) ''' and ''' inputname(2) ''' dynSeries object must have at least one common date!'])
end

if a.init<b.init
    n = b.init-a.init;
    b.data = [NaN(n,b.vobs); b.data];
    b.nobs = b.nobs+n;
    b.init = init;
end

if a.init>b.init
    n = a.init-b.init;
    a.data = [NaN(n,a.vobs); a.data];
    a.nobs = a.nobs+n;
    a.init = init;
end

if last_a>last_b
    n = last_a-last_b;
    b.data = [b.data; NaN(n,b.vobs)];
    b.nobs = b.nobs+n;
    return
end

if last_a<last_b
    n = last_b-last_a;
    a.data = [a.data; NaN(n,a.vobs)];
    a.nobs = a.nobs+n;
    return
end

%@test:1
%$ % Define a datasets.
%$ A = rand(8,3); B = rand(7,2);
%$
%$ % Define names
%$ A_name = {'A1';'A2';'A3'};
%$ B_name = {'B1';'B2'};
%$
%$ % Define initial dates
%$ A_init = '1990Q1';
%$ B_init = '1989Q2';
%$
%$ % Instantiate two dynSeries objects
%$ ts1 = dynSeries(A,A_init,A_name);
%$ ts2 = dynSeries(B,B_init,B_name);
%$
%$ try
%$   [ts1, ts2] = align(ts1, ts2);
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$ 
%$ if t(1)
%$   t(2) = dyn_assert(ts1.nobs,ts2.nobs);
%$   t(3) = dyn_assert(ts1.init==ts2.init,1);
%$   t(4) = dyn_assert(ts1.data,[NaN(3,3); A], 1e-15);
%$   t(5) = dyn_assert(ts2.data,[B; NaN(4,2)], 1e-15);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a datasets.
%$ A = rand(8,3); B = rand(7,2);
%$
%$ % Define names
%$ A_name = {'A1';'A2';'A3'};
%$ B_name = {'B1';'B2'};
%$
%$ % Define initial dates
%$ A_init = '1990Q1';
%$ B_init = '1990Q1';
%$
%$ % Instantiate two dynSeries objects
%$ ts1 = dynSeries(A,A_init,A_name);
%$ ts2 = dynSeries(B,B_init,B_name);
%$
%$ try
%$   [ts1, ts2] = align(ts1, ts2);
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dyn_assert(ts1.nobs,ts2.nobs);
%$   t(3) = dyn_assert(ts1.init==ts2.init,1);
%$   t(4) = dyn_assert(ts1.data,A, 1e-15);
%$   t(5) = dyn_assert(ts2.data,[B; NaN(1,2)], 1e-15);
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define a datasets.
%$ A = rand(8,3); B = rand(7,2);
%$
%$ % Define names
%$ A_name = {'A1';'A2';'A3'};
%$ B_name = {'B1';'B2'};
%$
%$ % Define initial dates
%$ A_init = '1990Q1';
%$ B_init = '1990Q1';
%$
%$ % Instantiate two dynSeries objects
%$ ts1 = dynSeries(A,A_init,A_name);
%$ ts2 = dynSeries(B,B_init,B_name);
%$
%$ try
%$   [ts2, ts1] = align(ts2, ts1);
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dyn_assert(ts1.nobs,ts2.nobs);
%$   t(3) = dyn_assert(ts1.init==ts2.init,1);
%$   t(4) = dyn_assert(ts1.data,A, 1e-15);
%$   t(5) = dyn_assert(ts2.data,[B; NaN(1,2)], 1e-15);
%$ end
%$ T = all(t);
%@eof:3