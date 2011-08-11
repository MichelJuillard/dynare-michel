function a = horzcat(b,c)

%@info:
%! @deftypefn {Function file} {@var{a} =} horzcat (@var{b},@var{c})
%! @anchor{horzcat}
%! @sp 1
%! Method of the dynSeries class.
%! @sp 1
%! Merge two Dynare time series class. This method overloads the horizontal concatenation operator, so that
%! two time series objects can be merged using the following syntax
%!
%!     a = [b, c];
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item b
%! Dynare time series object, instantiated by @ref{dynSeries}.
%! @item c
%! Dynare time series object, instantiated by @ref{dynSeries}.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @var
%! @item a
%! Dynare time series object.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @ref{descriptive_statistics}
%!
%! @strong{This function calls:}
%! @ref{dynSeries}
%!
%! @strong{Remark 1.} It is assumed that the two time series objects have the same frequencies. The two time series objects can cover
%! different time ranges.
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011 Dynare Team
% stephane DOT adjemian AT univ DASH lemans DOT fr
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

if ~isa(b,'dynSeries')
    error('dynSeries::horzcat: First input argument has to be a Dynare time series object!')
end

if ~isa(c,'dynSeries')
    error('dynSeries::horzcat: Second input argument has to be a Dynare time series object!')
end

if b.freq ~= c.freq
    error('dynSeries::horzcat: Two time series objects must have common frequency!')
else
    a = dynSeries();
    a.freq = b.freq;
end

d_nobs_flag = 0;
if b.nobs ~= c.nobs
    % error('dynSeries::horzcat: Two time series objects must have the same number of observations!')
    d_nobs_flag = 1;
else
    a.nobs = b.nobs;
end

d_init_flag = 0;
if isequal(b.init,c.init)
    a.init = b.init;
else
    % error('dynSeries:: Two time series objects must have common initial date!')
    % set a.init equal to min(b.init,c.init)
    if b.init(1)<c.init(1)
        d_init_flag = 1;
        a.init = b.init;
    elseif b.init(1)==c.init(1)
        if b.init(2)<c.init(2)
            d_init_flag = 1;
            a.init = b.init;
        else
            d_init_flag = 2;
            a.init = c.init;
        end
    else
        d_init_flag = 2;
        a.init = c.init;
    end
end

d_last_flag = 0;
if isequal(b.last,c.last)
    a.last = b.last;
else
    % error('dynSeries:: Two time series objects must have common initial date!')
    % set a.last equal to max(b.last,c.last)
    if b.last(1)<c.last(1)
        d_last_flag = 2;
        a.last = c.last;
    elseif b.last(1)==c.last(1)
        if b.last(2)<c.last(2)
            d_last_flag = 2;
            a.last = c.last;
        else
            d_last_flag = 1;
            a.last = b.last;
        end
    else
        d_last_flag = 1;
        a.last = b.last;
    end
end

a.vobs = b.vobs+c.vobs;
a.name = char(b.name,c.name);
a.tex  = char(b.tex,c.tex);

if ~( d_nobs_flag(1) || d_init_flag(1) || d_last_flag(1) )
    a.time = b.time;
    a.data = [b.data,c.data];
else
    [junk,ib] = setdiff(b.time,c.time,'rows');
    [junk,ic] = setdiff(c.time,b.time,'rows');
    [junk,jb,jc] = intersect(b.time,c.time,'rows');
    a.time = [b.time(ib,:); b.time(jb,:); c.time(ic,:)];
    a.time = sortrows(a.time,[1 2]);
    a.nobs = rows(a.time);
    a.data = NaN(a.nobs,a.vobs);
    [junk,ia,ib] = intersect(a.time,b.time,'rows');
    a.data(ia,1:b.vobs) = b.data;
    [junk,ia,ic] = intersect(a.time,c.time,'rows');
    a.data(ia,b.vobs+(1:c.vobs)) = c.data;
end


%@test:1
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$ B = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = char('A1','A2');
%$ B_name = char('B1','B2');
%$
%$ % Define expected results.
%$ e.time = [transpose(1:10), ones(10,1)];
%$ e.freq = 1;
%$ e.name = char('A1','A2','B1','B2');
%$ e.data = [A,B];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$ ts2 = dynSeries(B,[],B_name,[]);
%$
%$ % Call the tested method.
%$ ts3 = [ts1,ts2];
%$
%$ % Check the results.
%$ t(1) = dyn_assert(ts3.time,e.time);
%$ t(2) = dyn_assert(ts3.freq,e.freq);
%$ t(3) = dyn_assert(ts3.data,e.data);
%$ t(4) = dyn_assert(ts3.name,e.name);
%$ T = all(t);
%@eof:1

%@test:2
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$ B = [transpose(5:12),2*transpose(5:12)];
%$
%$ % Define names
%$ A_name = char('A1','A2');
%$ B_name = char('B1','B2');
%$
%$ % Define initial date
%$ A_init = 2001;
%$ B_init = 2005;
%$
%$ % Define expected results.
%$ e.time = [transpose(2000+(1:12)), ones(12,1)];
%$ e.freq = 1;
%$ e.name = char('A1','A2','B1','B2');
%$ e.data = [ [A; NaN(2,2)], [NaN(4,2); B]];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dynSeries(A,A_init,A_name,[]);
%$ ts2 = dynSeries(B,B_init,B_name,[]);
%$
%$ % Call the tested method.
%$ ts3 = [ts1,ts2];
%$
%$ % Check the results.
%$ t(1) = dyn_assert(ts3.time,e.time);
%$ t(2) = dyn_assert(ts3.freq,e.freq);
%$ t(3) = dyn_assert(ts3.data,e.data);
%$ t(4) = dyn_assert(ts3.name,e.name);
%$ T = all(t);
%@eof:2

%@test:3
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:7),2*transpose(1:7)];
%$ B = [transpose(5:11),2*transpose(5:11)];
%$
%$ % Define names
%$ A_name = char('A1','A2');
%$ B_name = char('B1','B2');
%$
%$ % Define initial date
%$ A_init = '1950Q1';
%$ B_init = '1950Q3';
%$
%$ % Define expected results.
%$ e.time = [ 1950, 1; 1950, 2; 1950, 3; 1950, 4; 1951, 1; 1951, 2; 1951, 3; 1951, 4; 1952, 1];
%$ e.freq = 4;
%$ e.name = char('A1','A2','B1','B2');
%$ e.data = [ [A; NaN(2,2)], [NaN(2,2); B]];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dynSeries(A,A_init,A_name,[]);
%$ ts2 = dynSeries(B,B_init,B_name,[]);
%$
%$ % Call the tested method.
%$ ts3 = [ts1,ts2];
%$
%$ % Check the results.
%$ t(1) = dyn_assert(ts3.time,e.time);
%$ t(2) = dyn_assert(ts3.freq,e.freq);
%$ t(3) = dyn_assert(ts3.data,e.data);
%$ t(4) = dyn_assert(ts3.name,e.name);
%$ T = all(t);
%@eof:3

%@test:4
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:7),2*transpose(1:7)];
%$ B = [transpose(5:9),2*transpose(5:9)];
%$
%$ % Define names
%$ A_name = char('A1','A2');
%$ B_name = char('B1','B2');
%$
%$ % Define initial date
%$ A_init = '1950Q1';
%$ B_init = '1950Q3';
%$
%$ % Define expected results.
%$ e.time = [ 1950, 1; 1950, 2; 1950, 3; 1950, 4; 1951, 1; 1951, 2; 1951, 3];
%$ e.freq = 4;
%$ e.name = char('A1','A2','B1','B2');
%$ e.data = [ A, [NaN(2,2); B]];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dynSeries(A,A_init,A_name,[]);
%$ ts2 = dynSeries(B,B_init,B_name,[]);
%$
%$ % Call the tested method.
%$ ts3 = [ts1,ts2];
%$
%$ % Check the results.
%$ t(1) = dyn_assert(ts3.time,e.time);
%$ t(2) = dyn_assert(ts3.freq,e.freq);
%$ t(3) = dyn_assert(ts3.data,e.data);
%$ t(4) = dyn_assert(ts3.name,e.name);
%$ T = all(t);
%@eof:4
